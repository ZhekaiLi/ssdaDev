% case 263
function [ ssda_results, probdata ] = Test_ss_da_GPinMatlab_concentratedCompare(lsf,probdata,analysisopt,gfundata,femodel,randomfield)


global nfun           % # LSF evaluations
global net gpopt ep;  % net: GP model; gpopt: optimizing parameter; 
                      % ep: expectation maximization options
global data;
global gpm;
global GPMs;

GPMs = {};

ssda_restart_from_step = analysisopt.ssda_restart_from_step; % -inf

if ssda_restart_from_step < 0
   nfun = 0;

   flag_plot = analysisopt.flag_plot;
   flag_plot_gen = analysisopt.flag_plot_gen;
   flag_cov_pf_bounds = analysisopt.flag_cov_pf_bounds;

   % Extract model data
   marg = probdata.marg;
   R = probdata.correlation;
   transf_type = probdata.transf_type; % 3

   % Find number of random variables
   nrv = size(marg,1);
  
   echo_flag = 1;

   block_size = analysisopt.block_size;

   rand_generator = analysisopt.rand_generator;
   stdv1 = analysisopt.stdv_sim;   % std for the crude MC 

   analysisopt.num_sim = 300; % intial 1000
   num_sim = analysisopt.num_sim;  % # samples used in each level
   
   ssda_Data.pf_target = analysisopt.pf_target; % 0.1
   ssda_Data.rho = analysisopt.rho; % 0.8
   ssda_Data.Pc = analysisopt.Pc; % 1

   point = zeros(nrv,1);

   % Modify correlation matrix and perform Cholesky decomposition
   if ~isfield(probdata,'Lo')
      % Compute corrected correlation coefficients
      [ Ro, dRo_drho, dRo_dthetafi, dRo_dthetafj ] = mod_corr_solve( marg, R , 0);
      probdata.Ro = Ro;

      % Cholesky decomposition
      % Lo = (chol(Ro))'; probdata.Lo = Lo;
      [Lo,ierr] = my_chol(Ro); if  ierr>0, return, end
      probdata.Lo = Lo;

      iLo = inv(Lo);
      probdata.iLo = iLo;
   end
   % Establish Cholesky decomposition of covariance matrix in crude MC
   chol_covariance = stdv1 * eye(nrv);  % chol_covariance = chol(covariance)
end

%% Subset Simulation - Step 0 (Monte Carlo Simulation)
ssda_Data.Nb_step = 0; %nb_step = 0;
ssda_Data.Seed = twister('state')';

% Initializations
k = 0;
percent_done = 0;

ssda_Data.U = zeros(nrv,num_sim);
ssda_Data.G = zeros(1,num_sim);

if ssda_restart_from_step < 0
   
   while k < num_sim

      block_size = min( block_size, num_sim-k );
      k = k + block_size;

      % Generate realizations of random U-vector
      allu = point*ones(1,block_size) + chol_covariance * inv_norm_cdf(twister(nrv,block_size));

      ssda_Data.U(:,(k-block_size+1):k) = allu;

      % Transform into original space
      allx = u_to_x(allu,probdata);
      ssda_Data.X(:,(k-block_size+1):k) = allx;

      % Evaluate limit-state function
      allg = gfun(lsf,allx,'no ',probdata,analysisopt,gfundata,femodel,randomfield);
      ssda_Data.G((k-block_size+1):k) = allg;
                 
      if floor( k/num_sim * 20 ) > percent_done
         percent_done = floor( k/num_sim * 20 );
         if echo_flag
            fprintf(1,'Subset step #%d - %d%% complete\n',ssda_Data.Nb_step,percent_done*5);
         end
      end

   end

   ssda_Data.Neval   = num_sim;       % number of LSFs evaluation
   ssda_Data.N       = num_sim;
   ssda_Data.Indices = [1 num_sim];
   ssda_Data.y       = [];            % threshold sequence
   ssda_Data.Indgerm = [];
   ssda_Data.AccRate = [];
   ssda_Data.p       = [];            % intermediate failure probability sequence
   if flag_cov_pf_bounds == 1
      ssda_Data.cov_pf_step  = [];
   end

   % Find y-threshold value.
   % Carry out intermediate calculations for the final estimation
   % of the coefficient variation of the failure probability pf.
   ssda_Data = ssda_y_threshold(ssda_Data,analysisopt);

   ssda_Data.Seed = [ ssda_Data.Seed twister('state')' ];
  
   % Train GP in matlab
   data.x = ssda_Data.U;
   data.y = ssda_Data.G;
   
   gpm = fitrgp(data.x', data.y');
   GPMs{1, 1} = gpm;
end

if ssda_Data.y ~= 0

   % Subset Simulation - Step 1
   ssda_Data.Nb_step = ssda_Data.Nb_step + 1; % nb_step = SubsetData.Nb_step;
   
   % Subset Simulation - Step 1 - Run simulations
   ssda_Data = Test_ssda_step_GPinMatlab_concentratedCompare(ssda_Data,num_sim,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
end

%% Loop until y-threshold equal to zero
while ssda_Data.y(end) ~= 0
   if ssda_restart_from_step < 0
      % Find y-threshold value.
      % Carry out intermediate calculations for the final estimation
      % of the coefficient variation of the failure probability pf.
      ssda_Data = ssda_y_threshold(ssda_Data,analysisopt);
   end
 
   ssda_Data.Seed = [ ssda_Data.Seed twister('state')' ];
   
   if ssda_Data.y(end) == 0, break; end
   
   % when the last three y is same, the model is not renewable
   if ssda_Data.y(end) == ssda_Data.y(end-1) && ssda_Data.y(end-1) == ssda_Data.y(end-2), break; end  

   % Subset Simulation - Step > 1
   ssda_Data.Nb_step = ssda_Data.Nb_step + 1; % nb_step = ssda_Data.Nb_step;

   ssda_Data = Test_ssda_step_GPinMatlab_concentratedCompare(ssda_Data,num_sim,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
end

% Assess approximation bounds for the coeficient of variation of the failure probability pf
if analysisopt.flag_cov_pf_bounds == 1
   ssda_Data.cov_pf_bounds = [ sqrt(sum((ssda_Data.cov_pf_step).^2)) sqrt(sum(sum((ssda_Data.cov_pf_step)'*(ssda_Data.cov_pf_step)))) ];
end

% Plot generalized reliability index vs. threshold value, with +/-2 stdv interval, for each subset step.
% Normal and lognormal hypothesis based on first subset step (CMC) samples
if echo_flag
   cov_pf_bounds = [];
   Pf = [];
   Beta = []; Beta_inf = []; Beta_sup = [];
   for nb_step = 0:ssda_Data.Nb_step
      cov_pf_inf = sqrt(sum((ssda_Data.cov_pf_step(1:(nb_step+1))).^2));
      cov_pf_sup = sqrt(sum(sum((ssda_Data.cov_pf_step(1:(nb_step+1)))'*(ssda_Data.cov_pf_step(1:(nb_step+1))))));
      cov_pf_bounds = [ cov_pf_bounds [ cov_pf_inf ; cov_pf_sup ] ];
      pf = prod(ssda_Data.p(1:(nb_step+1)));
      Pf = [ Pf pf ];
      Beta  = [ Beta -inv_norm_cdf(pf) ];
      Beta_inf = [ Beta_inf -inv_norm_cdf(pf+2*cov_pf_inf*pf) ];
      Beta_sup = [ Beta_sup -inv_norm_cdf(pf-2*cov_pf_inf*pf) ];
      if nb_step == 0
         mean0 = mean(ssda_Data.G(1:ssda_Data.Indices(1,2)));
         stdv0 = std(ssda_Data.G(1:ssda_Data.Indices(1,2)));
      	cov0 = stdv0/mean0;
      	zeta0 = (log(1+cov0^2))^0.5;
      	lambda0 = log(mean0) - 0.5*zeta0^2;
      end
   end
   [ Pf-2*cov_pf_bounds(1,:).*Pf ; Pf; Pf+2*cov_pf_bounds(1,:).*Pf ];
   [ Beta_sup; Beta; Beta_inf ];
   figure
   hbeta = errorbar(ssda_Data.y,Beta,Beta-Beta_inf,Beta_sup-Beta);
%    hold on
%    hnorm = plot(ssda_Data.y,-(ssda_Data.y-mean0)/stdv0,'r--');
%    hlogn = plot(ssda_Data.y,-(log(ssda_Data.y)-lambda0)/zeta0,'r:');
   set(hbeta,'LineWidth',2)
%    set(hnorm,'LineWidth',2)
%    set(hlogn,'LineWidth',2)
   set(gca,'FontSize',14);
   xlabel('{\ity_i}','FontSize',14);
   ylabel('{\it\beta_i}','FontSize',14);
end


ssda_Data.pf = prod(ssda_Data.p);

ssda_results.pf         = ssda_Data.pf;
ssda_results.beta       = -inv_norm_cdf(ssda_Data.pf);
ssda_results.nfun       = nfun;
ssda_results.ssda_Data  = ssda_Data;
ssda_results.net        = net;
ssda_results.gpopt      = gpopt;

if analysisopt.flag_cov_pf_bounds == 1
   ssda_results.cov_pf  = ssda_Data.cov_pf_bounds;
end
