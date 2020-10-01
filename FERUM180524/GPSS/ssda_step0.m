function ssda_Data = ssda_step(ssda_Data,num_sim,lsf,probdata,analysisopt,gfundata,femodel,randomfield)

% Perform a single subset step

global nfun

nrv = size(probdata.marg,1);  % dim(u)

rand_generator = analysisopt.rand_generator;
if isfield(analysisopt,'echo_flag')
   echo_flag = analysisopt.echo_flag;
else
   echo_flag = 1;
end
flag_plot      = analysisopt.flag_plot;
flag_plot_gen  = analysisopt.flag_plot_gen;
 
rho    = ssda_Data.rho;  % correlation coefficient in GCS
Indgerm  = ssda_Data.Indgerm(end,1:length(find(ssda_Data.Indgerm(end,:))));
subgermU = ssda_Data.U(:,Indgerm);  % seeds for intermediate sampling
subgermG = ssda_Data.G(Indgerm);    % LSF value for seeds
Nseeds = size(subgermG,2);

subsetU   = [];
subsetG  = [];

Nb_generation = 0;
MHrate = [];        % acceptance rate

while size(subsetU,2) < num_sim
  
   Nb_generation = Nb_generation + 1;

   subtempu = ss_gcs(subgermU,rho,rand_generator); % Gaussian conditional sampling

   % Initializations
   block_size = analysisopt.block_size;
   
   k = 0;
   percent_done = 0;

%    Nseeds = I_eval;
%    subtempu = subtempu;
%    subtempg = subgermG;
%    subtempg = subtempg(1,I_eval);
      
   while k < Nseeds
   
      block_size = min( block_size, Nseeds-k );
      k = k + block_size;
      
      allu = subtempu(:,(k-block_size+1):k);

      % Transform into original space
      allx = u_to_x(allu,probdata);

      % Evaluate limit-state function
      allG = gfun(lsf,allx,'no ',probdata,analysisopt,gfundata,femodel,randomfield);
      
      subtempg((k-block_size+1):k) = allG;
      
      % Train GP

      if floor( k/Nseeds * 20 ) > percent_done
         percent_done = floor( k/Nseeds * 20 );
        if echo_flag
            fprintf(1,'Subset step #%d - Generation #%d - %d%% complete\n',ssda_Data.Nb_step,Nb_generation,percent_done*5);
        end
      end
      
   end
   
   
   ssda_Data.Neval = ssda_Data.Neval + Nseeds;            % # LSFs evaluation
   ssda_Data.N     = ssda_Data.N     + size(subtempu,2);
   
   % Selection of points not in the failure domain
   ind = find( subtempg > ssda_Data.y(end) );  % rejected points

   if nrv == 2 && flag_plot == 1 && flag_plot_gen == 1
      ssda_Data.Usubtemp = subtempu;
      ssda_Data.ind = ind;
      ssda_Data.Nb_generation = Nb_generation;
      ssda_Data = ssda_graph(ssda_Data,2);
   end
   
   MHrate = [ MHrate ; 100*(Nseeds - length(ind))/Nseeds ];

   % replace the rejected points with original points
   subtempu(:,ind) = subgermU(:,ind);
   subtempg(ind) = subgermG(ind);
   % collect the generated samples
   subsetU = [ subsetU subtempu ]; 
   subsetG = [ subsetG subtempg ];
   % Update the seeds
   subgermG = subtempg;
   subgermU = subtempu;
  
end

ssda_Data.U       = [ ssda_Data.U subsetU ];
ssda_Data.G       = [ ssda_Data.G subsetG ];
ssda_Data.Indices = [ ssda_Data.Indices; ssda_Data.Indices(end,2)+1 ssda_Data.Indices(end,2)+size(subsetU,2) ];
ssda_Data.AccRate = [ ssda_Data.AccRate [ Nb_generation; mean(MHrate) ] ];
