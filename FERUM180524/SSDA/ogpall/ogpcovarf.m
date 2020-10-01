function covf = ogpcovarf(x1, x2,varargin)
%OGPCOVARF Calculate the covariance function for the OGP.
%
%	Description
%
%	COVF = OGPCOVARF(X1, X2) considers the global OGP data structure NET
%	together with two matrices X1, X2 and computes the matrix of
%	covariance function COVF.  If called with a single input X1, the
%	function returns only the diagonal of OGPCOVARF(X1,X1).
%
%	If the output is not one-dimensional, then the returned matrix has
%	NOUT times the lenght of the inputs X1 and X2 respectively.  Using a
%	single input set gives the covariance function as DIM(X1),NOUT,NOUT
%	three dimensional matrix.
%
%	The scalar product (or distances) between inputs X1 and X2 are
%	computed using the scaling from EXP(NET.KPAR(1:NET.NIN)). For all
%	covariance types there is an amplitude EXP(NET.KPAR(NET.NIN+1)).
%	Other kernel-specific parameters can be specified at higher
%	positions: NET.KPAR(NET.NIN+2) and above.  The number of parameters
%	in the kernel function can be arbitrarily high.
%
%	The available kernel functions are the following:
%
%	 'SQEXP'  - squared exponential.
%
%	 'POLYEXP'  - polynomial exponential kernel. A product of the
%	exponential and the polynomial kernels.
%
%	 'RATQUAD'  - rational quadtratic kernel. The order of the kernel
%	function is NET.KPAR(2).
%
%	 'POLY'  - polynomial kernel.  The order of the polynomial is
%	EXP(NET.KPAR(2)).
%
%	 'LINSPLINE'  - linear spline kernel.     'SINC'  - the sinc kernel
%	with frequency cutoff at   EXP(NET.KPAR(2)).
%
%	 'FSIN'  - the FSIN kernel used in density estimation.
%
%	 'USER'  - user-specified kernel function.
%
%	If the user-specified covariance function is used, i.e.
%	NET.COVARFN='USER', NET.KPAR still contains the hyperparameters - the
%	parameters of the kernel, but there are two additional fields in the
%	GLOBAL structure NET:
%
%	 KFNADDR  - the address of the function that   returns the covariance
%	kernel.
%
%	 GRADKADDR  - the address of the function returning the gradient   of
%	the kernel w.r.to the kernel parameters.
%
%
%	An example is COV_MATERN which implements the Matern covariance
%	function.
%
%	See the online (HTML) documentation for details
%
%	See also
%	OGP, OGPINIT, OGPCOVARP, OGPCOVDIAG, COV_MATERN
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% consistency checking
errstring = consist(net, 'ogp', x1);
if ~isempty(errstring);
  error(errstring);
end

%% ------- returning the diagonal of the kernel matrix --------
if nargin==1;		      % returns the diagonal of
  covf = ogpcovdiag(x1);
  return;
end;

%  -----------  returning a matrix  ----------

if size(x1, 2) ~= size(x2, 2)
  error('Number of variables in x1 and x2 must be the same');
end

beta = spdiags(exp(net.kpar(1:net.nin)'),...
	       0, ...
	       net.nin*net.nout,net.nin*net.nout);
aA   = exp(net.kpar(1+net.nin));
ord  = exp(net.kpar(2+net.nin));

n1 = size(x1, 1);
n2 = size(x2, 1);
covf = zeros(n1*net.nout,n2*net.nout);
if ~n1 | ~n2; return; end;

switch upper(net.covarfn)
 case 'SQEXP'		      % Squared exponential
  % Compute the weighted squared distances between x1 and x2
  covf = (x1.*x1)*beta*ones(net.nin, n2) - 2*x1*beta*x2' ... 
	 + ones(n1, net.nin)*beta*(x2.*x2)';
  covf = full(aA * exp( - 0.5*covf));
  
 case 'OU'		      % exponential Ornstein-Uhlenbeck
  for iD=1:net.nin;
    covf = covf + ...
	   sqrt(beta(iD,iD))*abs(x1(:,iD)*ones(1,n2) - ones(n1,1)*x2(:,iD)');
  end;
  covf = aA*exp(-covf);
  
 case 'POLYEXP'
  % Compute the weighted squared distances between x1 and x2
  covf = (x1.*x1)*beta*ones(net.nin, n2) - 2*x1*beta*x2' ... 
	 + ones(n1, net.nin)*beta*(x2.*x2)';
  covf = sqrt(covf);
  covf = aA*exp(-covf.^2./2).*(1+covf).^ord;
  
 case 'RATQUAD'		      % Rational quadratic
  beta = beta./ord;
  % Compute the weighted squared distances between x1 and x2
  covf = (x1.*x1)*beta*ones(net.nin, n2) - 2*x1*beta*x2' ... 
	 + ones(n1, net.nin)*beta*(x2.*x2)';

  covf = aA * ((1 + covf).^(-ord));

 case 'POLY'		      % Polynomial
  covf    = 1 + x1*beta*x2';
  
  covf = aA *covf.^round(ord);
  
 case 'LINSPLINE'	      % linear splines
  % K(x,y) = min(x,y)^3/3 + min(x,y)^2*abs(x-y)/2 + x*y + 1   
  % THE MULTIDIM. case is the sum of the individual kernels.
  if find([x1(:) x2(:)] < 0);
    error('Nodes for LINSPLINE need to be positive');
  else
    tPar = exp(net.kpar((net.nin+1):end));
    if length(tPar==1);
      tPar(2:net.nin) = tPar;
    end;
    for iIn=1:net.nin;
      matX1 = diag(x1(:,iIn)) * ones(n1,n2);
      matX2 = ones(n1,n2) * diag(x2(:,iIn));
      minXX = min(matX1,matX2);
      tCov  = beta(iIn,iIn)^(3/2)*(minXX.^3/3+minXX.^2.*abs(matX1-matX2)/2) ...
	      + beta(iIn,iIn)*matX1.*matX2 + 1;
      covf  = covf + tPar(iIn) * tCov;
    end;
  end;

 case 'SINC';		      % SINC kernel for density estimation
  tPar = exp(net.kpar((net.nin+1):end));
  if length(tPar==1);
    tPar(2:net.nin) = tPar;
  end;
  for iIn=1:net.nin;
    [t1,t2] = meshgrid(x2(:,iIn),x1(:,iIn));
    covf = covf + tPar(iIn) * sinc(beta(iIn,iIn)*(t2-t1));
  end;

 case 'FSIN';		      % FINITE Dimensional periodic kernel
  for iIn = 1:net.nin;
    [t1,t2] = meshgrid(x2(:,iIn),x1(:,iIn));
    tCov    = beta(iIn,iIn)*(t2 - t1);
    % constraining the inputs to [0,1]
    tCov  = tCov - floor(tCov);
    U     = aA * (2*ord+1)*ones(size(tCov));
    nZ    = find(tCov);
    U(nZ) = aA * sin(pi*(2*ord+1)*tCov(nZ))./sin(pi*tCov(nZ));
    covf  = covf + U;
  end;

 case 'USER';		      % user-specified covariance function
  covf = feval(net.kfnaddr,x1,x2);
  
 otherwise
  error(['Unknown covariance function ', net.covarfn]);
end;
