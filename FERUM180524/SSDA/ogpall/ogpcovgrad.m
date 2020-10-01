function g = ogpcovgrad(x1,covf,iPar)
%OGPCOVGRAD Evaluate gradient for kernel parameters (except bias) for the Sparse GP
%
%	Description
%
%	G = OGPCOVGRAD(NET) takes a Sparse OGP data structure NET and
%	evaluates the gradient G of the negative log-likelihood with respect
%	to the hyperparameters of the model. The output G is a matrix with
%	columns of length #BV^2 (the kernel matrix is put into vector
%	format).  Each column in G corresponds to a hyper-parameter and the
%	order is provided by OGPPAK.
%
%	G = OGPCOVGRAD(NET,XTRAIN) additionally to the data structure NET,
%	the procedure takes the elements of the training set XTRAIN and
%	returns on the columns the derivatives of the kernel function.  The
%	result is thus a matrix with (NBV X NTR) lines and the number of
%	columns the number of hyperparameters.
%
%	If the kernel is specified by the user - i.e. kernel type is 'USER' -
%	then the last group of parameters is given by the user and there must
%	exist a field 'FNGRAD' in the structure 'NET.KPAR'.  For details on
%	specifying a different kernel see the documentation for function
%	OGPCOVARF.
%
%	See also
%	OGP, OGPPAK, OGPEVIDGRAD, OGPTRAIN
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% To speed up computations, we might give covf - the kernel matrix.

% !!! COVF must be a matrix of size [n^2 1]
% !!! COVF must have the bias term subtracted !!!

% checking parameters
if nargin<1;  x1   = []; end;
if nargin<2;  covf = []; end;
if nargin<3;  iPar = 0;
else
  if length(iPar)>1;
    fprintf(['Single gradient computation implemented, you asked:\n['...
	     '%d]'],iPar);
    error;
  end
end;

% gradient for the USER-defined covariance matrix
if strcmp(upper(net.covarfn),'USER');
  g = feval(net.gradkaddr,x1,covf,iPar);
  return;
end;

% IF X1 is empty, BV is assigned to it
if ~prod(size(x1));  x1 = net.BV; end;

% computing the kernel matrix
if ~prod(size(covf));	      % COVF was empty
  covf = ogpcovarf(net.BV,x1);
end;

% setting temporal variables
nx    = size(net.BV,1);
n1    = size(x1,1);
nAll  = nx * n1;

% preparing the inputs
beta = spdiags(...
    reshape(exp(net.kpar(1:net.nin)), [net.nin 1]), ...
    0,net.nin,net.nin);
aA   = exp(net.kpar(1+net.nin));
ord  = exp(net.kpar(2+net.nin));

% !! SAVING MEMORY !!
if strcmp(upper(net.covarfn),'POLY')
  % covR - the reduced kernel....
  covR = zeros(nAll,1);
  ord   = round(ord);
  if ~mod(ord,2)	      % EVEN numbers
    for iDim=1:net.nin;
      xt = repmat(x1(:,iDim),[nx 1]);
      % FIRST:  [1:n,1:n,...,1:n]
      zt = reshape(repmat(net.BV(:,iDim)', [n1 1]),...
		   [1 nAll])';
      covR = covR + full((xt.*zt).*beta(iDim,iDim));
    end;
    covR = aA*(1+covR).^(ord-1);
  else
    covR = exp(log(abs(covf))*((ord-1)./ord)+net.kpar(1+net.nin)/ord);
  end;
end;

nin   = net.nin; 	      %    ZERO == ALL
iLS   = [1:nin];
% we need the inputs only for the gradients w.r.t length-scale
if ~iPar | (iPar<=net.nin);
  if iPar;		      % NONZERO == SINGLE
    nin = 1;        iLS = iPar;
  end;
  % PUT: [1,1,1, ... , ... , n,n,n]
  x1   = repmat(x1(:,iLS),[nx 1]);
  % FIRST:  [1:n,1:n,...,1:n]
  z    = reshape(repmat(net.BV(:,iLS)', [n1 1]),...
		 [nin nAll])';
end;

% computing gradients
switch upper(net.covarfn)
 case 'SQEXP'		      % Squared exponential
  if ~iPar | iPar<=net.nin;
    z    = full((x1 - z).^2 * beta(iLS,iLS)./2);
  end;
  % ALL COLUMNS RETURNED
  if ~iPar;		      % DEFAULT: !all!
    % the order of building 'g' is defined by ogppak
    % first the input scales
    g(:,1:nin)  = - repmat(covf,[1 nin]).* z;
    % then the amplitude 
    g(:,net.nin +1) = covf;
    g(:,net.nin +2) = zeros(nAll,1);
  else			      % SELECTED COLUMNS ONLY
    if iPar<=net.nin;
      g = - covf .* z;
    elseif iPar==net.nin+1;
      g = covf;
    else
      g = zeros(nAll,1);
    end;
  end;

 case 'OU'  		      % ORNSTEIN-UHLENBECK
  if ~iPar | iPar<=net.nin;
    beta = spdiags(...
	reshape(exp(net.kpar(1:net.nin)/2), [net.nin 1]), ...
	0,net.nin,net.nin);
    z    = full(abs(x1 - z) * beta(iLS,iLS) );
  end;

  % OPTIONALLY THERE IS A SINGLE COMPONENT RETURNED
  if ~iPar;		      % DEFAULT: !all!
    % the order of building 'g' is defined by ogppak
    % first the scaling parameters
    g(:,1:net.nin)  = - repmat(covf,[1 net.nin]).* z/2;
    % then the amplitude
    g(:,net.nin +1) = covf;
    g(:,net.nin +2) = zeros(nAll,1);
  else
    if iPar<=net.nin;
      g = - covf .* z ./ 2;
    elseif iPar==net.nin+1;
      g = covf;
    else
      g = zeros(nAll,1);
    end;
  end;

 case 'RATQUAD'		      % Rational quadratic
  if ~iPar | iPar<=net.nin;
    % Compute the weighted squared distances
    z     = full((x1 - z).^2 * beta(iLS,iLS));
  end
  % OPTIONALLY THERE IS A SINGLE COMPONENT RETURNED
  if ~iPar;		      % DEFAULT: !all!
    log1sZ= -log(covf/aA)./ord;
    % the order of building 'g' is defined by ogppak
    % first the scaling parameters
    g(:,1:net.nin)  = - repmat(covf.*exp(-log1sZ),[1 net.nin]).*z;
    % (amplitude)
    g(:,net.nin+1)  = covf;
    % (order)
    g(:,net.nin+2)  = - covf .* ord .* ...
	( log1sZ - 1 + exp(-log1sZ));
  else
    if iPar==net.nin+1;
      g = covf;
      return;
    end;
    log1sZ= -log(covf./aA)./ord;
    if iPar<=net.nin;
      g = - covf.*exp(-log1sZ).*z;
    elseif iPar==net.nin+2;
      g = - covf .* ord .* (log1sZ - 1 + exp(-log1sZ));
    end;
  end;
    
 case 'POLY'		      % Polynomial
  if ~iPar | iPar<=net.nin;
    z      = full((x1 .* z) * beta(iLS,iLS));
  end;
  % OPTIONALLY THERE IS A SINGLE COMPONENT RETURNED
  if ~iPar;		      % DEFAULT: !all!
    % the order of building 'g' is defined by ogppak/ogpunpak
    % first the scaling parameters
    g(:,1:net.nin)  = repmat(covR,[1 net.nin]) .* z * ord;
    % (amplitude)
    g(:,net.nin+1)  = covf;
    % (order)
    g(:,net.nin+2)  = zeros(nAll,1);
  else
    if iPar<=net.nin;
      g = covR .* z * ord;
    elseif iPar==net.nin+1;
      g = covf;
    else
      g = zeros(nAll,1);
    end;
  end;

 otherwise
    error(['Gradient not implemented for the covariance function ', ...
	   net.covarfn]);
end;
