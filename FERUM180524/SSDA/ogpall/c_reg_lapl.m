function [logLik,K1,K2,varargout] = c_reg_lapl(likpar,y,mX,sigX2,varargin)
%C_REG_LAPL Online update coeffs. for regression with Laplace noise.
%
%	Description
%
%	[LOGLIK,K1,K2] = C_REG_LAPL(LIKPAR,Y,MX,VARX) - returns the
%	coefficients for the online update of the GP regression with Gaussian
%	noise assumption.
%
%	Parameters:
%
%	 LIKPAR  - likelihood parameters (stored in NET).
%
%	 Y  - observed noisy output.
%
%	 MX  - mean of the GP marginal at the new input.
%
%	 VARX  - the variance of the GP marginal at the new input.
%
%	 LOGLIK  - the logarithm of the averaged likelihood function.
%
%	 K1,K2  - the first and second derivatives of LOGLIK.
%
%	The noise parameter LAMBDA=LIKPAR is set at the initialisation.
%
%	See also
%	OGP, OGPTRAIN, OGPFWD, ERR_ABS
%

%	Copyright (c) Lehel Csato (2001-2004)

varargout = {[]};

threshold = 1e-10;
% checking the validity of the variance for the GP marginal
if sigX2 < eps;
  fprintf('Warning! Pred. distr. var.: %f\n',sigX2);
  sigX2 = eps;
end;

sK1   = - sign(y-mX);
sqSig = sqrt(sigX2);
sS    = sqSig/likpar;   % \lambda*sigma
y     = abs((y - mX)./sqSig); % the function is symmetric

% calling the auxiliary function that perfomrs the approximations.
l_n = slapint(sS,y);
l_d = slapint(sS,-y);

if l_n-l_d<-3;
  logLik = exp(l_n-l_d);
else
  logLik = log(1+exp(l_n-l_d));
end;
logLik = -log(likpar*2) + l_d + logLik;

if nargout==1; return; end;   % STOP if nothing else is wanted

% compute derivatives
rat = exp(l_n - l_d);

K1     = sK1.*(1 - 2./(1+rat))./likpar;

K2     = - 2*exp(-y.^2/2 - l_d)./sqrt(2*pi)./sqSig./(1+rat)./likpar + ...
	 4./(1+rat).*(1-1./(1+rat))./likpar.^2; 

% adjusting K2 - to make the INFERENCE alg. more stable
tolL  = 1e6;
indBad = find(-K2<1./tolL);
if length(indBad);
  K2(indBad) = - 1./tolL;
end;

%%%%%%%%%%%%% INTERNAL FUNCTION - SLAPINT - INTERNAL FUNCTION %%%%%%%%%%
function [l_n]  = slapint(s,x)
% function [l_n, l_d] = slapint(s,x) - is an auxiliary function required
% for the posterior average over a laplace likelihood.

% the log-average has two terms: g(s,x) and g(s,-x). This function
% computes  log(g(s,x)).  The high-order (8) asymptotics are NEEDED
% since for large values of s=\lambda*sigma_x the two values are close
% to zero - and we need to obtain information about the second
% derivative.

% where the approximations should start
thresh = 6;

% selecting where to do the asymptotic approx.
% iAs -- the indices for the inputs which are APPROXIMATED
iAs = find(-x-s<-thresh);
% iN -- the indices for inputs that use normal computation
iN  = setdiff([1:length(x)],iAs);

iS = 1; iNS = 1;
if length(s)>1;
  iS = iAs; iNS = iN;
end;
% result
l_n = zeros(size(x));

tt       = 1./(x(iAs)+s(iS)).^2;
l_n(iAs) = - x(iAs).^2./2 - log((x(iAs)+s(iS))*sqrt(2*pi)) - ...
	(1 - (5/2 - (37./3 - 353/4*tt).*tt).*tt).*tt;
tt       = x(iN);
l_n(iN)  = s(iN).*(s(iN)/2+tt) + log((1-erf((tt+s(iN))./sqrt(2)))./2);
