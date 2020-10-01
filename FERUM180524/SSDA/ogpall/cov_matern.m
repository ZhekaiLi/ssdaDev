function covf = cov_matern(x1,x2);
%COV_MATERN returns the Matern covariance function.
%
%	Description
%
%	COVF = COV_MATERN(X1, X2) takes two matrices X1, X2 of input vectors
%	and computes the matrix of covariance function COVF.  If called with
%	a single input X1, the function returns the diagonal of
%	COV_MATERN(X1,X1).
%
%	This kernel function assumes that the output is one-dimensional.
%
%	See the online (HTML) documentation for details
%
%	See also
%	OGPCOVARF, OGPCOVGRAD, COVGRAD_MATERN, DEMOGP_MATERN
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% argument checking
if nargin < 1;
  error('Wrong number of parameters in COV_MATERN');
end;

% consistency checking
errstring = consist(net, 'ogp', x1);
if ~isempty(errstring);
  error(errstring);
end

[n1, dim] = size(x1);
sigma = exp(net.kpar(net.nin+1));
nu    = exp(net.kpar(net.nin+2));
beta  = exp(net.kpar(1:net.nin)');

if nargin==1;		      % diagonal elements only
  covf = sigma * repmat(eye(net.nout),[n1 1]);
else  			      % computing a full matrix
  % consistency checking
  errstring = consist(net, 'ogp', x2);
  if ~isempty(errstring);
    error(errstring);
  end
  n2   = size(x2,1);
  covf = zeros(n1,n2);
  for iD = 1:dim;
    covf = covf + ...
	   (repmat(x1(:,iD),1,n2)-repmat(x2(:,iD)',n1,1)).^2*beta(iD);
  end;
  covf = sqrt(abs(full(covf)));
  covf = matern(nu,sigma,covf);
end;

