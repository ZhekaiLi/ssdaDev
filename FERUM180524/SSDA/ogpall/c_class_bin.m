function [loglik,K1,K2,varargout] = c_class_bin(likpar,y,x,sigX2,varargin)
%C_CLASS_BIN Online update coeffs. for binary classification.
%
%	Description
%
%	[LOGLIK,K1,K2] = C_CLASS_BIN(LIKPAR,Y,MX,VARX,PM) - returns the
%	coefficients for the online update of the GP assuming classification.
%
%	Parameters:
%
%	 LIKPAR  - the likelihood parameter (output noise).
%
%	 Y  - observed noisy labels.
%
%	 MX  - mean of the GP marginal at the new input.
%
%	 VARX  - the variance of the GP marginal at the new input.
%
%	 LOGLIK  - the logarithm of the averaged likelihood function.
%
%	 K1,K2  - the first and second derivatives of LOGLIK.
%
%	Although in practice not important, the noise variance is
%	net.likpar(1) which should be set to a nonzero value when
%	initialising the model.
%
%	See also
%	OGP, OGPTRAIN, ERR_2CLASS, DEMOGP_CLASS
%

%	Copyright (c) Lehel Csato (2001-2004)

varargout = {[]};

% checking the validity of the variance for the GP marginal
if sigX2 < eps;
  fprintf('Warning! Pred. distr. var.: %f\n',sigX2);
  sigX2 = eps;
end;

sigX2 = sigX2 + likpar(1);
x = y * x ./ sqrt(sigX2);

% avoiding numerical errors: splitting the inputs...
iSm    = find(x<-5);
iNn    = find(x>=-5 & x<6);
iLl    = find(x>=6);

tERF   = (1 + erf(x(iNn)/sqrt(2)))/2;
tt     = x(iSm).^(-2);

if nargout>1;		      % compute only if required

  tExpX  = exp(-x(iNn).^2/2)/sqrt(2*pi);
  c      = tExpX./tERF;

  K1(iSm) = - y(iSm).*x(iSm).*(1 + (1-2*tt).*tt)./ sqrt(sigX2(iSm));
  K2(iSm) = (-1 + (1-6*tt).*tt) ./ sigX2(iSm);

  K1(iLl) = 0;
  K2(iLl) = 0;

  K1(iNn) = y(iNn) .* c ./ sqrt(sigX2(iNn));
  K2(iNn)  = -c.*(x(iNn) + c) ./ sigX2(iNn);
end;

% only -LOGPRED
loglik(iSm) = -.5*(log(2*pi) + 1./tt) - log(-x(iSm)) - ...
    (1-(2.5-(37/3-353/4.*tt).*tt).*tt).*tt;
loglik(iLl) = 0;
loglik(iNn) = log(tERF);
