function [ ssda_results, probdata ] = ss_da0(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

%     Finite Element Reliability Using Matlab, FERUM, Version 4.0, 2009 
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     A copy of the GNU General Public License is found in the file
%     <gpl-3.0.txt> following this collection of FERUM program files.
%     This license can be also found at: <http://www.gnu.org/licenses/>.    
%
%     For more information on FERUM, visit: <http://www.ifma.fr/FERUM/>

global nfun           % # LSF evaluations
global net gpopt ep;  % net: GP model; gpopt: optimizing parameter; 
                      % ep: expectation maximization options
gpini(analysisopt);   % initialize the GP model

if isfield(analysisopt,'ssda_restart_from_step')
   ssda_restart_from_step = analysisopt.ssda_restart_from_step;
else
   ssda_restart_from_step = -inf;
end

%% Load data from previous subset step, if necessary ( ss_restart_from_step >= 0 )
if ssda_restart_from_step >= 0
   eval(['load ssda_restart_from_step_' num2str(ssda_restart_from_step) '.mat']);
   twister('state',ssda_Data.Seed(:,ssda_restart_from_step+2)');
end


if ssda_restart_from_step < 0

   nfun = 0;

   if ~isfield(analysisopt,'flag_plot')
      analysisopt.flag_plot = 0;
   end
   flag_plot = analysisopt.flag_plot;
   if ~isfield(analysisopt,'flag_plot_gen')
      analysisopt.flag_plot_gen = 0;
   end
   flag_plot_gen = analysisopt.flag_plot_gen;
   if ~isfield(analysisopt,'flag_cov_pf_bounds')
      analysisopt.flag_cov_pf_bounds = 1;
   end
   flag_cov_pf_bounds = analysisopt.flag_cov_pf_bounds;

   % Extract model data
   marg = probdata.marg;
   R = probdata.correlation;
   transf_type = probdata.transf_type;

   % Find number of random variables
   nrv = size(marg,1);
   
   if isfield(analysisopt,'echo_flag')
      echo_flag = analysisopt.echo_flag;
   else
      echo_flag = 1;
   end

   block_size = analysisopt.block_size;

   rand_generator = analysisopt.rand_generator;
   stdv1 = analysisopt.stdv_sim;   % std for the crude MC 
   num_sim = analysisopt.num_sim;  % # samples used in each level
   
   if ~isfield(analysisopt,'pf_target')
       analysisopt.pf_target = 0.1;
   end
   ssda_Data.pf_target = analysisopt.pf_target;
   
   if ~isfield(analysisopt,'rho')
      analysisopt.rho = 0.8;
   end
   ssda_Data.rho = analysisopt.rho;
   
   if ~isfield(analysisopt,'Pc')
       analysisopt.Pc = 1;
   end
   ssda_Data.Pc = analysisopt.Pc;

   point = zeros(nrv,1);

   % Modify correlation matrix and perform Cholesky decomposition
   if ~isfield(probdata,'Lo')

      if transf_type == 3

         % Compute corrected correlation coefficients
         switch probdata.Ro_method
            case 0
               Ro = mod_corr( marg, R );
            case 1
               [ Ro, dRo_drho, dRo_dthetafi, dRo_dthetafj ] = mod_corr_solve( marg, R , 0);
         end
         probdata.Ro = Ro;

         % Cholesky decomposition
         % Lo = (chol(Ro))'; probdata.Lo = Lo;
         [Lo,ierr] = my_chol(Ro); if  ierr>0, return, end
         probdata.Lo = Lo;

         iLo = inv(Lo);
         probdata.iLo = iLo;

      end

   end

   % Establish Cholesky decomposition of covariance matrix in crude MC
   chol_covariance = stdv1 * eye(nrv);  % chol_covariance = chol(covariance);

end


%% Subset Simulation - Step 0 (Monte Carlo Simulation)

ssda_Data.Nb_step = 0; %nb_step = 0;
ssda_Data.Seed = twister('state')';

if ssda_restart_from_step == -1
   save ssda_restart_from_step_0.mat ssda_Data
end

if ssda_restart_from_step < 0

   % Initializations
   k = 0;
   percent_done = 0;

   ssda_Data.U = zeros(nrv,num_sim);
   ssda_Data.G = zeros(1,num_sim);
   
end

if nrv == 2 && exist('flag_plot','var') == 1 && flag_plot == 1
   close all
end

if ssda_restart_from_step < 0
   
   while k < num_sim

      block_size = min( block_size, num_sim-k );
      k = k + block_size;

      % Generate realizations of random U-vector
      switch rand_generator
         case 0  % Matlab function
            allu = point*ones(1,block_size) + chol_covariance * randn(nrv,block_size);
          otherwise  % twister
            allu = point*ones(1,block_size) + chol_covariance * inv_norm_cdf(twister(nrv,block_size));
      end

      ssda_Data.U(:,(k-block_size+1):k) = allu;

      % Transform into original space
      allx = u_to_x(allu,probdata);

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
   
   if ssda_restart_from_step == -1
      ssda_restart_from_step0 = ssda_restart_from_step;
      analysisopt0 = analysisopt;
      clear ss_restart_from_step
      analysisopt = rmfield(analysisopt,'ss_restart_from_step');
      save ss_restart_from_step_0.mat
      ssda_restart_from_step = ssda_restart_from_step0;
      analysisopt = analysisopt0;
   end
   
end

%% Train GP
Ns = 800;
U = lhsnorm(zeros(nrv,1),eye(nrv),Ns,'off');
% Transform into original space
X = u_to_x(U',probdata);
% Evaluate limit-state function
G = gfun(lsf,X,'no ',probdata,analysisopt,gfundata,femodel,randomfield);
traingp(U,G');
%    ytest = ogpfwd(ssda_Data.U');

%% 
% plot
if nrv == 2 && exist('flag_plot','var') == 1 && flag_plot == 1
   ssda_Data = ssda_graph(ssda_Data,0);
end


if ssda_Data.y ~= 0

   % Subset Simulation - Step 1
   
   ssda_Data.Nb_step = ssda_Data.Nb_step + 1; % nb_step = SubsetData.Nb_step;
   
   if ssda_Data.Nb_step > ssda_restart_from_step
      if ~( ssda_restart_from_step < 0 )
         analysisopt.ss_restart_from_step = -1;
         ssda_restart_from_step = -1;
      end
   end

   if ssda_restart_from_step < 0
      % Subset Simulation - Step 1 - Run simulations
      ssda_Data = ssda_step(ssda_Data,num_sim,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
   end

   if nrv == 2 && exist('flag_plot','var') == 1 && flag_plot == 1
      ssda_Data = ssda_graph(ssda_Data,1);
   end

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

   if ssda_restart_from_step == -1
      ssda_restart_from_step0 = ssda_restart_from_step;
      analysisopt0 = analysisopt;
      clear ssda_restart_from_step
      analysisopt = rmfield(analysisopt,'ssda_restart_from_step');
      eval(['save ssda_restart_from_step_' num2str(ssda_Data.Nb_step) '.mat']);
      ssda_restart_from_step = ssda_restart_from_step0;
      analysisopt = analysisopt0;
   end
   
   if nrv == 2 && exist('flag_plot','var') == 1 && flag_plot == 1
      ssda_Data = ssda_graph(ssda_Data,0);
   end
   
   if ssda_Data.y(end) == 0, break; end

   % Subset Simulation - Step > 1
   ssda_Data.Nb_step = ssda_Data.Nb_step + 1; % nb_step = ssda_Data.Nb_step;
   
   if ssda_Data.Nb_step > ssda_restart_from_step
      if ~( ssda_restart_from_step < 0 )
         analysisopt.ssda_restart_from_step = -1;
         ssda_restart_from_step = -1;
      end
   end

   if ssda_restart_from_step < 0
      % Subset Simulation - Step > 1 - Run simulations
      ssda_Data = ssda_step(ssda_Data,num_sim,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
   end
   
   if nrv == 2 && exist('flag_plot','var') == 1 && flag_plot == 1 && ssda_Data.y(end) ~= 0
      ssda_Data = ssda_graph(ssda_Data,1);
   end
  
end

% Assess approximation bounds for the coeficient of variation of the failure probability pf
if analysisopt.flag_cov_pf_bounds == 1
   ssda_Data.cov_pf_bounds = [ sqrt(sum((ssda_Data.cov_pf_step).^2)) sqrt(sum(sum((ssda_Data.cov_pf_step)'*(ssda_Data.cov_pf_step)))) ];
end

if nrv == 2 && exist('flag_plot','var') == 1 && flag_plot == 1
   ssda_Data = ssda_graph(ssda_Data,3);
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
   hold on
   hnorm = plot(ssda_Data.y,-(ssda_Data.y-mean0)/stdv0,'r--');
   hlogn = plot(ssda_Data.y,-(log(ssda_Data.y)-lambda0)/zeta0,'r:');
   set(hbeta,'LineWidth',2)
   set(hnorm,'LineWidth',2)
   set(hlogn,'LineWidth',2)
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
