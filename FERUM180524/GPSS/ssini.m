
global nfun           % # LSF evaluations
% global net gpopt ep;  % net: GP model; gpopt: optimizing parameter; 
                      % ep: expectation maximization options

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
