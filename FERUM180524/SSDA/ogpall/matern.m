function covf = matern(nu,sigma,ddd);
%MATERN	Computes the Matern kernel
%
%	Description
%
%	K = MATERN(NU, SIGMA, D) computes the value of the Matern kernel of
%	order NU and scales the output with factor SIGMA.
%
%	The inputs NU and SIGMA have to be scalars and the output has the
%	dimensionality of the third input argument.  The implementation uses
%	the scaled version of BESSELK which has a more stable behaviour for
%	large orders of the Bessel function.
%
%	See also
%	COV_MATERN, COVGRAD_MATERN, DEMOGP_MATERN
%

%	Copyright (c) Lehel Csato (2001-2004)

% checking parameters
error(nargchk(3, 3, nargin));
if prod(size(nu))~=1 | prod(size(sigma))~=1;
  error('Only scalar paratemers are allowed in MATERN');
end;
if nu<0;
  error('Order of Bessel function has to be positive');
end;

if nu<100;		      % compute the Bessel kernel
  
  w      = warning;
  warning('off');
  ddd    = sqrt(2*nu)*ddd;
  logT   = log(ddd);
  logCov = log(besselk(nu, ddd, 1)) - ddd + nu*logT;
  logCov = logCov + log(sigma)-(nu-1)*log(2)-gammaln(nu);
  covf   = exp(logCov);
  covf(isnan(covf)) = sigma;
  covf(isinf(covf)) = sigma;
  % For large degree NU, BESSELK will give Inf. Z==0 may cancel this
  besselInf = find(isinf(logCov));
  if any(ddd(besselInf)>0),
    warning(sprintf(['Some of the Bessel functions are Infinite. Degree\n'...
		     'of Matern function NU is too large: %f'],nu));
  end;
  warning(w);

else			      % SUBSTITUTE the approximation
  covf = sigma.*exp(-ddd.^2/2);
end;
