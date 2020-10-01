%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'PROBDATA'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc
% Names of random variables. Default names are 'x1', 'x2', ..., if not explicitely defined.
probdata.name =  { 'x1'
                   'x2'};

% Marginal distributions for each random variable
% probdata.marg =  [ (type) (mean) (stdv) (startpoint) (p1) (p2) (p3) (p4) (input_type); ... ];								  
probdata.marg =  [  1  0  1  0  nan  nan  nan  nan 0 ;
                    1  0  1  0  nan  nan  nan  nan 0 ];

% Correlation matrix
probdata.correlation = eye(2);
                                                                                               

probdata.transf_type = 3;
probdata.Ro_method   = 1;
% Flag for computation of sensitivities w.r.t. thetag parameters of the limit-state function
% 1: all sensitivities assessed, 0: no sensitivities assessment
probdata.flag_sens   = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'ANALYSISOPT'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

analysisopt.multi_proc           = 1;        % 1: block_size g-calls sent simultaneously
                                             %    - gfunbasic.m is used and a vectorized version of gfundata.expression is available.
                                             %      The number of g-calls sent simultaneously (block_size) depends on the memory
                                             %      available on the computer running FERUM.
                                             %    - gfunxxx.m user-specific g-function is used and able to handle block_size computations
                                             %      sent simultaneously, on a cluster of PCs or any other multiprocessor computer platform.
                                             % 0: g-calls sent sequentially
analysisopt.block_size           = 10000;  % Number of g-calls to be sent simultaneously


% Subset Simulation w/ delayed acceptance (SS_DA) analysis options
analysisopt.num_sim              = 1000;     % Number of samples (MC,IS), number of samples per subset step (SS) or number of directions (DS)
analysisopt.rand_generator       = 1;        % 0: default rand matlab function, 1: Mersenne Twister (to be preferred)
analysisopt.sim_point            = 'origin'; % 'dspt': design point, 'origin': origin in standard normal space (simulation analysis)
analysisopt.stdv_sim             = 1;        % Standard deviation of sampling distribution in simulation analysis

analysisopt.rho                  = 0.8;      % Correlation coefficient
analysisopt.pf_target            = 0.2;      % Target probability for each subset step
analysisopt.flag_cov_pf_bounds   = 1;        % 1: calculate upper and lower bounds of the coefficient of variation of pf, 0: no calculation 
analysisopt.ssda_restart_from_step = -inf;   % i>=0 : restart from step i, -inf : all steps, no record (default), -1 : all steps, record all
analysisopt.Pc                   = 1;        % accepting probability
analysisopt.flag_plot            = 1;        % 1: plots at each step (2 r.v. examples only), 0: no plots
analysisopt.flag_plot_gen        = 1;        % 1: intermediate plots for each MCMC chain (2 r.v. examples only), 0: no plots

% GP surrogate
% GPmodel.ktype                = 'sqexp';  % kernal type: sqexp, ratquad, poly, matern
% GPmodel.dim                  = 2;        % dimension of RVs
% GPmodel.Nc                   = 1e3;      % size of the maximal training set
% GPmodel.epic                 = 1e-3;     % innovation threshold
% analysisopt.GPmodel          = GPmodel;  % GP surrogate setting

% Subset Simulation (SS) analysis options
analysisopt.width                = 2;        % Width of the proposal uniform pdfs
analysisopt.pf_target            = 0.2;      % Target probability for each subset step
analysisopt.flag_cov_pf_bounds   = 1;        % 1: calculate upper and lower bounds of the coefficient of variation of pf, 0: no calculation 
analysisopt.ss_restart_from_step = -inf;     % i>=0 : restart from step i, -inf : all steps, no record (default), -1 : all steps, record all
analysisopt.flag_plot            = 0;        % 1: plots at each step (2 r.v. examples only), 0: no plots
analysisopt.flag_plot_gen        = 0;        % 1: intermediate plots for each MCMC chain (2 r.v. examples only), 0: no plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'GFUNDATA' (one structure per gfun)  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Type of limit-state function evaluator:
% 'basic': the limit-state function is defined by means of an analytical expression or a Matlab m-function,
%          using gfundata(lsf).expression. The function gfun.m calls gfunbasic.m, which evaluates gfundata(lsf).expression.
% 'xxx':   the limit-state function evaluation requires a call to an external code.  The function gfun.m calls gfunxxx.m,
%          which evaluates gfundata(lsf).expression where gext variable is a result of the external code.
gfundata(1).evaluator  = 'basic';
gfundata(1).type       = 'expression';   % Do no change this field!

% Expression of the limit-state function:
gfundata(1).expression = 'b - x2 - kappa * ( x1 - e ).^2';

% cg deterministic parameters of the limit-state function (optional, may also be defined in probdata.marg)
% gfundata(1).cg = [ (cg1) (cg2) ... ]
%               or [ (cg1) (cg2) ... ]'
gfundata(1).cg         = [ 5 0.5 0.1 ];

% Names of cg parameters (cg1, cg2, ... if not defined):
% gfundata(1).cgname = {'name1' 'name2' ... }
%                   or {'name1' 'name2' ... }'
gfundata(1).cgname     = { 'b' 'kappa' 'e' };

% Flag for computation of sensitivities w.r.t. thetag parameters of the limit-state function
% 1: all sensitivities assessed, 0: no sensitivities assessment
gfundata(1).flag_sens  = 0;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'FEMODEL'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

femodel = [];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  DATA FIELDS IN 'RANDOMFIELD'  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

randomfield = [];