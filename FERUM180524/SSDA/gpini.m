function gpini(analysisopt,probdata,gfundata,femodel,randomfield)
% Initialization of the GP surrogate

GP = analysisopt.GPmodel;

global net gpopt ep;

% INITIALISING THE GP
if strcmp(GP.ktype,'sqexp')
    ogp(GP.dim,1,'sqexp',...          % input, output dimension, kernel type
    log([ones(1,GP.dim)*0.1 10 1]));   % length-scales, amp. & order
elseif strcmp(GP.ktype,'ratquad')
    ogp(GP.dim,1,'ratquad',...        % input, output dimension, kernel type
    log([ones(1,GP.dim)*0.1 10 2]));   % length-scales, amp. & order
elseif strcmp(GP.ktype,'poly')
    ogp(GP.dim,1,'poly',...           % input, output dimension, kernel type
    log([ones(1,GP.dim)*0.1 40 4]));   % length-scales, amp. & order
elseif strcmp(GP.ktype,'matern')
    ogp(GP.dim,1,'user',...           % input, output dimension, kernel type
    log([0.04 0.04 0.15 0.15 0.01 0.01 0.01 0.002 40 5/2])); % length-scales, amp. & order
%     log([ones(1,GP.dim)*0.1 10 3])); % length-scales, amp. & order
    net.kfnaddr   = @cov_matern;      % KERNEL
    net.gradkaddr = @covgrad_matern;  % GRADIENT wrt parameters
end

% HYPERPRIORS used in hyperparameter optimisation.
ogphypcovpar(5e-2);                   % hyplambda specifies the precision (inverse variance) of the parameters

% ASSIGNING A LIKELIHOOD FUNCTION
ogpinit(@c_reg_gauss,1e-4);    	      % Gaussian likelihood for regression
net.thresh    = GP.epic;              % the admission threshold for new BVs
net.maxBV     = GP.Nc;                % capacity of training set
net.isBVfixed = 1;

% Initialising GPOPT - optimisation and display options
gpopt = defoptions;      % NO log-average

% the control of HYPPAR optimisation
    gpopt.covopt.opt    = foptions;% default options
    gpopt.covopt.opt(1) = 0;  % display values
    gpopt.covopt.opt(9) = 0;  % CHECK/NOT gradients
    gpopt.covopt.opt(2) =1e-3;% precision in computing f. value
    gpopt.covopt.opt(3) =1e-3;% display values
    gpopt.covopt.opt(14)=6;% number of iterations
    gpopt.covopt.fnopt ='conjgrad';
%     
% gpopt.covopt.opt    = foptions;   % default options
% gpopt.covopt.opt     = 2;
% gpopt.covopt.fnopt   = 'scg';
gpopt.covopt.hyopt   = 1;         % = 1: optimize hyperparameter, = 0: use default hyperparameter
gpopt.postopt.isep   = 1;         % USING the TAP/EP algorithm
gpopt.postopt.itn    = 2;         % number of EP iteration with CHANGING BVs
gpopt.postopt.fixitn = 3;         % FIXING the BV set.
gpopt.postopt.isrm   = 1;         % flag to remove BV or not
gpopt.postopt.nem    = 3;         % number of em steps in optimizing the hyperparameters


%% Train GP
if ~isfield(GP,'net')
Ns = 300;
U = lhsnorm(zeros(GP.dim,1),eye(GP.dim),Ns,'off');
% Transform into original space
X = u_to_x(U',probdata);
% Evaluate limit-state function
G = gfun(1,X,'no ',probdata,analysisopt,gfundata,femodel,randomfield);
traingp(U,G');
else
    net = GP.net;
end
%    ytest = ogpfwd(ssda_Data.U');

