function ogptrain(xTrain,yTrain)
%OGPTRAIN Inference for Sparse Gaussian Processes.
%
%	Description
%
%	The function OGPTRAIN(X_TRAIN,Y_TRAIN) takes the GLOBAL Gaussian
%	Process data structure NET and performs the following operations:
%
%	 -  approximates the posterior process.
%
%	 -  adjusts covariance function parameters.
%
%	 -  adjusts likelihood function parameters.
%
%
%	The above steps constitute a single cycle in the EM algorithm built
%	for the joint optimisation of the posterior process and the
%	hyperparameters.
%
%	Notice that - due to the changes in the hyperparameters - the
%	posterior process at the end of OGPTRAIN is no longer optimal.  Thus,
%	if prediction is wanted, then one should perform an extra computation
%	step using OGPPOST.  This is not done in OGPTRAIN.
%
%	The calculations are governed by the fields of the GLOBAL structure
%	GPOPT: the indicators influencing calculation of the posterior are
%	grouped into the sub-structure GPOPT.POSTOPT.
%
%	Given an approximation to the posterior process (step 1), we then
%	adjust the parameters of the covariance kernel and the likelihood.
%	The covariance parameters are optimised using a conjugate gradient
%	algorithm, the specific algorithm and the number of steps can be
%	altered via the structure GPOPT.COVOPT.
%
%	The optimisation of the likelihood parameters is done as the last
%	step of the EM procedure.  It is independent of the optimisation of
%	covariance function parameters.  It can be done using gradients (like
%	G_L_GAUSS) or EM (see EM_GAUSS).
%
%	See also
%	OGP, OGPPOST, OGPPAK, OGPEVID, OGPEVIDGRAD
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% THIS FUNCTION IS A WRAPPER 
% 
% which includes:
%  1. computing the approximation to the posterior with fixed hyp.par.
%  2. computing a new set of covariance parameters based on the post.
%  3. computing the new likelihood parameters

%%%%%%%%%%%%%%%%%%% ARGUMENT CHECKING %%%%%%%%%%%%%%%%%%%%%%%%%%
nTrain    = size(xTrain,1);
% handling models without direct output
if strcmpi(net.outtype,'DIRECT'); 
  errstring = consist(net, 'ogp', xTrain, yTrain);
elseif strcmpi(net.outtype,'INDIRECT');
  errstring = consist(net, 'ogp', xTrain);
  if size(xTrain,1) ~=size(yTrain,1);
    errstring = ['Number of input patterns ', num2str(size(xTrain,1)), ...
		 ' does not match number of output patterns ', ...
		 num2str(size(yTrain,1))];
    error(errstring);
  end;
else
  errstring = [];
end;
if ~isempty(errstring);
  error(errstring);
end;
if nargin ~= 2;
  error('Insufficient number of arguments in calling OGPTRAIN');
end;
if isempty(gpopt);
  gpopt = defoptions;
end;

if isnumeric(gpopt);
  itn          = gpopt;
  gpopt     = defoptions;
  gpopt.postopt.itn = itn;
end;

if isfield(gpopt,'ptest') && gpopt.ptest;	
  % there has to be a test set!
  if ~(isfield(gpopt,'xtest') && isfield(gpopt,'ytest'));
    error('Wrong GPOPT: xtest or ytest not found');
  else 
    errstring = consist(net,'ogp',gpopt.xtest,gpopt.ytest);
    if ~isempty(errstring);
      error(errstring);
    end;
  end;
  if ~isfield(gpopt,'erraddr'); % default error function: quadratic
    gpopt.erraddr = @err_mse;
  end;
  % asking for test error if there are no observables is irreal!
  if ~strcmpi(net.outtype,'DIRECT');
    error(['When no direct observations, cannot compute test error -\n' ...
	   '     WRONG GPOPT ARGUMENT']);
  end;
end;
if ~isfield(gpopt,'freq');
  gpopt.freq = 1;
end;

%%%%%%%%%%%%%%%  APPROXIMATING THE POST. PROCESS  %%%%%%%%%%%%%%

ogppost(xTrain,yTrain);

% EXTRA STEPS performed to STABILISE the EP parameters
if gpopt.postopt.fixitn;
  oldInd        = net.isBVfixed;
  oldItn        = gpopt.postopt.itn;
  net.isBVfixed = 1;
  
  gpopt.postopt.itn = gpopt.postopt.fixitn;
  ogppost(xTrain,yTrain);

  net.isBVfixed     = oldInd;
  gpopt.postopt.itn = oldItn;
end;

%%%%%%%%%%%%%%%  POST==PRIOR   =>  SHORTEN THE LENGTHSCALE
if isempty(net.KB);
    net.kpar(1:net.nin) = net.kpar(1:net.nin) - 0.22314355;
    net.kpar(1+net.nin) = net.kpar(1+net.nin) + 0.22314355;
    net.likpar          = net.likpar/2;
%%%%%%%%%%%%%%% LIKELIHOOD PARAMETER OPTIMISATION %%%%%%%%%%%%%%
elseif ~isempty(net.likoptfn);% if there is an optimiser
  [postM, postV] = ogpfwd(xTrain);
  % adjusting likelihood parameters
  net.likpar = feval(net.likoptfn,net.likpar,...
		     yTrain,ep.cavM,ep.cavV,...
		     postM,postV);
end;

%%%%%%%%%%%%%%% COVARIANCE PARAMETER OPTIMISATION %%%%%%%%%%%%%%
% optimising COVARIANCE parameters
if sum(abs(net.w))<net.thresh;% if it makes sense ...
  return;
end;

% parameters for the (NETLAB) optimisation algorithm
netOpt      = foptions;
netOpt(1)   =   -1;	      % DISPLAY/NOT training values
netOpt(9)   =    0;	      % CHECK/NOT   gradients
netOpt(2)   = 1e-6;	      % precision in computing optimum values
netOpt(3)   = 1e-6;

% find the initial parameters
w0 = ogppak;

% setting the number of SCG iterations
if isempty(gpopt.covopt.opt); % no step size was given
  netOpt(14) = floor(2*length(w0));
elseif length(gpopt.covopt.opt) == 1;% only the stepsize was given
  netOpt(14) = gpopt.covopt.opt;
else			      % all other parameters were given;
  netOpt     = gpopt.covopt.opt;
end;

% apply optimisation algorithm
if ~isempty(gpopt.covopt.fnopt)
  [wf, netOpt] = feval(gpopt.covopt.fnopt, ...
		       'ogpevid',w0,netOpt,'ogpevidgrad');

  ogpunpak(wf);
  
  if netOpt(9);
    % gradchek(wf,'ogpevid','ogpevidgrad');
  end;
end;
