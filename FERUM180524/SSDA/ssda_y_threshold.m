function ssda_Data = ssda_y_threshold(ssda_Data,analysisopt)

% Evaluate the threshold value y of the intermediate limit-state-function
% Make preliminary calculations for estimating the coefficient of variation of the failure probability pf

p = ssda_Data.pf_target;

% define range for quantile analysis
IndInf = ssda_Data.Indices(end,1);
IndSup = ssda_Data.Indices(end,2);

% find values of LSF
G = ssda_Data.G(IndInf:IndSup);
N = length(G);

% y = ferum_quantile(G,p);  % get p-quantile of G
y = quantile(G,p);

if y < 0
   y = 0;
   p = sum( G < y ) / N;   % failure prob. in last step
end

ssda_Data.y = [ ssda_Data.y y ];  % threshold
ssda_Data.p = [ ssda_Data.p p ];  % failure probability

Indgerm = find( G < y );        % seeds for next simulation
% in case with repeated quantile
dif = analysisopt.num_sim*ssda_Data.pf_target - length(Indgerm);
if dif > 0;
    I_add = find( G == y );
    Indgerm = sort([Indgerm I_add(1:dif)]);
end
Indgerm = Indgerm + IndInf - 1; % index in the whole simulation data

% if y > 0
    if ssda_Data.Nb_step == 0
        ssda_Data.Indgerm = Indgerm;
    elseif length(Indgerm) > size(ssda_Data.Indgerm,2),
        ssda_Data.Indgerm = [ ssda_Data.Indgerm...
            zeros(size(ssda_Data.Indgerm,1),length(Indgerm)-size(ssda_Data.Indgerm,2)) ];
        ssda_Data.Indgerm = [ ssda_Data.Indgerm; Indgerm ];
    else
        ssda_Data.Indgerm = [ ssda_Data.Indgerm; Indgerm ];
    end
% end

% c.o.v. computation
if analysisopt.flag_cov_pf_bounds == 1
   
   if ssda_Data.Nb_step == 0

      cov_pf_step = sqrt( (1-p) / (p*N) );  % c.o.v for crude MC

   else   % c.o.v. for intermediate level

      Nc = length(find(ssda_Data.Indgerm(end-1,:)));
      N = IndSup - IndInf + 1;
      
      Rzero = p*(1-p);
      
      I =  G <= y;
      Indices = reshape((1:N),[],N/Nc)';
      gamma = 0;
      for k = 1:(N/Nc-1)
         Z = 0;
         for j = 1:Nc
            for l = 1:(N/Nc-k)
               Z = Z + I(Indices(l,j)) * I(Indices(l+k,j));
            end
         end
         rho(k) = ( 1/(N-k*Nc) * Z - p^2 ) / Rzero;
         gamma = gamma + 2 * (1-k*Nc/N) * rho(k);
      end

      cov_pf_step = sqrt( (1-p) / (p*N) * (1+gamma) );

   end
   
   ssda_Data.cov_pf_step = [ ssda_Data.cov_pf_step cov_pf_step ];  % c.o.v for each iteration
   
end