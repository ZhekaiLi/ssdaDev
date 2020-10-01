function ssda_Data = ssda_step_compareGPMs(ssda_Data,num_sim,lsf,probdata,analysisopt,gfundata,femodel,randomfield)

% Perform a single subset step

global nfun
global net gpopt ep;
global NETS;

nrv = analysisopt.GPmodel.dim;  % # RVs

if isfield(analysisopt,'echo_flag')
   echo_flag = analysisopt.echo_flag;
else
   echo_flag = 1;
end
flag_plot      = analysisopt.flag_plot;
flag_plot_gen  = analysisopt.flag_plot_gen;

rho    = ssda_Data.rho;             % correlation coefficient in GCS
Indgerm  = ssda_Data.Indgerm(end,1:length(find(ssda_Data.Indgerm(end,:))));
subgermU = ssda_Data.U(:,Indgerm);  % seeds for intermediate sampling
subgermG = ssda_Data.G(Indgerm);    % LSF value for seeds
Nseeds = size(subgermG,2);

subsetU     = subgermU;     % store generated U
subsetG     = subgermG;     % store generated G
net.Utrain  = [];     % store newly generated U
net.Gtrain  = [];     % store newly generated G

Nb_generation = 0; % iteration order for sequential sampling
GPrate = [];       % acceptance rate for GP surrogate
DArate = [];       % acceptance rate for the dalayed acceptance

% generate uniform random variable for comparison
rand_generator = analysisopt.rand_generator;
switch rand_generator
case 0
u01 = rand(Nseeds,num_sim/Nseeds-1);
otherwise
u01 = twister(Nseeds,num_sim/Nseeds-1);
end

Pr_old = ones(Nseeds,1);             % predicted failure probability in seeds

while size(subsetU,2) < num_sim
  
   Nb_generation = Nb_generation + 1;
   
   % First step of SS-DA
   % Gaussian conditional sampling
   subtempu = ss_gcs(subgermU,rho,rand_generator); 

   % Initializations
   block_size = analysisopt.block_size;
   
   k = 0;
   percent_done = 0;
  
   subtempg = zeros(1,Nseeds);          % store G
   Pr_fail = zeros(Nseeds,1);           % predicted failure probability
   I_reject = [];                       % indices for rejected samples
   I1_reject = [];                      % indices for rejected samples by GP
   Neval = 0;                           % # LSF evaulations
   % parallel computing
   while k < Nseeds
   
      block_size = min(block_size, Nseeds - k );
      k = k + block_size;
      subind = (k-block_size+1):k;
      
      % Second step of SS-DA
      allu = subtempu(:,subind);
      % Transform into original space
      allx = u_to_x(allu,probdata);
      
      % prediction by GP surrogate
      % ZK: 01/10
      % store the GP model in every step
      % compare the value calculate by different GP model that have been calculate
      % then choose the result with smallest variance
      meanG = [];
      varG = [];
      for i = 1:size(NETS, 2)
         net = NETS{1, i};
         [meanGtemp, varGtemp] = ogpfwd(allx');
         varGtemp(varG<0) = 0;
         if i == 1
            meanG = meanGtemp;
            varG = varGtemp;
         else
            better = find(varGtemp < varG);
            meanG(better) = meanGtemp(better);
            varG(better) = varGtemp(better);
         end
      end
      subtempg(subind) = meanG;  % store G
      % predicted failure probability
      Pr_fail(subind) = normcdf((ssda_Data.y(end)-meanG)./sqrt(varG));
      
      % acceptance ratio of GP
      a2 = Pr_fail(subind)./Pr_old(subind);
      % rejected samples by GP
      I2_reject = subind(a2<u01(subind,Nb_generation));
      % cumulative rejected samples
      I1_reject = [I1_reject I2_reject];
      
      % Third step of SS-DA         
      % sample indices needed to evaluate
      I_evalu = setdiff(subind(Pr_fail(subind)<=analysisopt.Pc),I2_reject);
      I_evalx = I_evalu - k + block_size;
      Neval = Neval + length(I_evalu);  % # LSF evaluations
      I3_reject = [];      % indices for samples rejected by LSF
      if ~isempty(I_evalu)
          % Evaluate limit-state function
          newG = gfun(lsf,allx(:,I_evalx),'no ',probdata,analysisopt,gfundata,femodel,randomfield);
          subtempg(I_evalu) = newG;  % store G
          % indices of rejected samples
          I3_reject = I_evalu(newG > ssda_Data.y(end));
          Pr_fail(setdiff(I_evalu,I3_reject)) = 1;  % update failure probability
          ogppost(allu(:,I_evalx)',newG');  % update GP
          % collect newly generated samples for GP training
          net.Utrain = [net.Utrain allu(:,I_evalx)];
          net.Gtrain = [net.Gtrain newG];
      end
      % cumulative rejected samples
      I_reject = sort([I_reject [I2_reject I3_reject]]);
      
      if floor( k/Nseeds * 20 ) > percent_done
         percent_done = floor( k/Nseeds * 20 );
        if echo_flag
            fprintf(1,'Subset step #%d - Generation #%d - %d%% complete\n',ssda_Data.Nb_step,Nb_generation,percent_done*5);
        end
      end
      
   end
   % end parallel computing
   
   ssda_Data.Neval = ssda_Data.Neval + Neval;            % # LSFs evaluation
   ssda_Data.N     = ssda_Data.N     + size(subtempu,2); % # samples
   
   % Selection of points not in the failure domain
   if nrv == 2 && flag_plot == 1 && flag_plot_gen == 1
      ssda_Data.Usubtemp = subtempu;
      ssda_Data.ind = I_reject;
      ssda_Data.Nb_generation = Nb_generation;
      ssda_Data = ssda_graph(ssda_Data,2);
   end
   
   % Acceptance rate by GP
   GPrate = [ GPrate ; (Nseeds - length(I1_reject))/Nseeds ];
   % Acceptance rate by SS-DA
   DArate = [ DArate ; (Nseeds - length(I_reject))/Nseeds ];

   % replace the rejected points with original points
   subtempu(:,I_reject) = subgermU(:,I_reject);
   subtempg(I_reject) = subgermG(I_reject);
   Pr_fail(I_reject) = Pr_old(I_reject);
   % collect the generated samples
   subsetU = [ subsetU subtempu ]; 
   subsetG = [ subsetG subtempg ];
   % Update the seeds
   subgermG = subtempg;
   subgermU = subtempu;
   Pr_old = Pr_fail;
  
end

% return results
ssda_Data.U       = [ ssda_Data.U subsetU ];
ssda_Data.G       = [ ssda_Data.G subsetG ];
ssda_Data.Indices = [ ssda_Data.Indices; ssda_Data.Indices(end,2)+1 ssda_Data.Indices(end,2)+size(subsetU,2) ];
ssda_Data.AccRate = [ ssda_Data.AccRate [ Nb_generation ; mean(GPrate); mean(DArate) ] ];

% train GP
if ~isempty(net.Utrain) && gpopt.covopt.hyopt == 0
%     ogpreset;
    ogptrain(net.Utrain',net.Gtrain');
    ogppost(net.Utrain',net.Gtrain');
    NETS{1, size(NETS, 2)+1} = net;
end
