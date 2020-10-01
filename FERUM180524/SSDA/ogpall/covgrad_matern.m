function g = covgrad_matern(x1,covf,iPar);
%COVGRAD_MATERN Returns the gradient of the Matern kernel.
%
%	Description
%
%	G = COVGRAD_MATERN computes the derivative of the Matern kernel.  If
%	an input XTRAIN is present, then it is assumed that the gradient of
%	the matrix COV_MATERN(NET.XTRAIN,NET.BV) is computed.
%
%	The kernel parameters are taken from the global variable NET.
%
%	With two parameters it is assumed that the covariance matrix of the
%	BV set and XTRAIN has been computed and a kernel matrix computation
%	is spared.
%
%	A further reduction in memory requirements is when we ask only a
%	single column of the gradient matrix, specified by a third argument:
%	g = covgrad_matern(xTrain,covf,iPar)
%
%	If IPAR==0 then all gradient elements are computed and returned.
%	Often this operation fails due to memory requirements.
%
%	There is no analytical formula to differentiate the Bessel function
%	with respect to its order, thus numerical (finite difference)
%	approximation is employed.
%
%	See also
%	OGPCOVGRAD, COV_MATERN, DEMOGP_MATERN
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% To speed up computations, we might give covf - the kernel matrix.

% !!! COVF must be a matrix of size [n^2 1]
% !!! COVF must have the bias term subtracted !!!

% checking parameters
if nargin<1;  x1   = []; end;
if nargin<2;  covf = []; end;
if nargin<3;  iPar = 0;  end;

% IF X1 is empty, then the BV is the second set
if ~prod(size(x1));
  x1 = net.BV;
end;

% setting temporal variables
nx    = size(net.BV,1);
n1    = size(x1,1);
nAll  = nx * n1;

% computing the parameters
beta   = reshape(exp(net.kpar(1:net.nin)), [net.nin 1]);
sigma  = exp(net.kpar(net.nin+1));
nu     = exp(net.kpar(net.nin+2));

% !! SAVING MEMORY !!
ddd    = zeros(nAll,1);
for iDim=1:net.nin;
  xt   = repmat(x1(:,iDim),[nx 1]);
  zt   = reshape(repmat(net.BV(:,iDim)',[n1 1]),[1 nAll])';
  if ~iPar;
    z(:,iDim) = (xt - zt).^2.*beta(iDim).*2.*nu;;
  elseif iDim == iPar;
    z  = (xt - zt).^2.*beta(iDim).*2.*nu;;
  end;
  ddd  = ddd + (xt - zt).^2.*beta(iDim).*2.*nu;
end;
ddd    = sqrt(abs(ddd));

% computing the kernel matrix
if prod(size(covf)) ~= prod(size(z)); % recompute COVF  
  covf = matern(nu,sigma,ddd./sqrt(2*nu));
end;

warning('off');

% computing derivatives w.r.t AMPLITUDE
if ~iPar;
  g(:,net.nin+1)  = covf;
elseif iPar==net.nin+1;
  g = covf;
  return;
end;

% the ratio of the Bessel functions
% needed in the gradients w.r.t length-scale
if nu<100;
  iS = find(ddd<1e-5);
  iL = setdiff([1:length(ddd)],iS);
  % computing for large distances
  dPlus(iL,:) = - 2*nu./ddd(iL) ...
      + besselk(nu+1,ddd(iL),1)./besselk(nu, ddd(iL), 1);
  dPlus(iS,:) = 0;
else
  dPlus = ddd./2./nu;
end;

% computing derivatives w.r.t input scaling
if iPar<=net.nin;
  rRatio                = ddd;
  rRatio(find(~rRatio)) = 1;
  rRatio                = z./repmat(rRatio,[1 net.nin])./2;
  if iPar>0;		      % single gradient is needed
    g = -covf.*dPlus.*rRatio(:,iPar);
    return;
  else
    g(:,1:net.nin)  = repmat(-covf.*dPlus,[1 net.nin]).*rRatio;
  end;
end;

% gradient w.r.t KERNEL order
if nu<100;
  tTerm = nu*(log(ddd/2) - psi(nu) ...
	      + dlogbessk(nu,ddd) ...
	      - ddd.*dPlus./2./nu);
  
  % we know that LIM(log(x/2) - psi(nu) + dlogbesselk(nu,x),x=+0)=0
  tTerm(find(abs(ddd)<1e-5)) = 0;
  if ~iPar;
    g(:,net.nin+2) =  covf.*(tTerm + 100000*eps);
  else			      % IT CAN ONLY ASK FOR ORDER
    g =  covf.*(tTerm + 100000*eps);
  end;
else
  if ~iPar;
    g(:,net.nin+2) = zeros(size(covf)) - 1e-10;
  else			      % IT CAN ONLY ASK FOR ORDER
    g = zeros(size(covf)) - 1e-10;
  end;
end;
warning('on')

%%%%%%%%%%%%%%%%%%%%%%%% ---   BUILTIN   --- %%%%%%%%%%%%%%%%%%%%%%%%%
function grBess = dlogbessk(nu, d)
% DLOGBESSK - Numerical gradient of LOG(BESSELK) with respect to degree NU
%
%   D = DLOGBESSK(NU, Z) Compute the derivate of the LOG- modified Bessel
%   function (second kind) with respect to the degree NU, evaluated at
%   Z. Since there is no closed form expression for this gradient, this is
%   computed numerically.
%   DLOGBESSK(NU, Z, DNU) uses a step size of DNU to find the numerical
%   gradient. Default value: NU*1E-6

% 
% Copyright (c) by Anton Schwaighofer (2003)
% $Revision: 1.1 $ $Date: 2003/01/27 19:40:50 $
% mailto:anton.schwaighofer@gmx.net
% 

error(nargchk(2,2,nargin));
dnu = nu*1e-4;

grBess = ( log(besselk(nu+dnu,d,1)) - log(besselk(nu-dnu,d,1)) ) ...
    ./dnu * 0.5;
