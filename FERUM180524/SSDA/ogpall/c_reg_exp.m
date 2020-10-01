function [logLik,K1,K2,varargout] = c_reg_exp(likpar,y,mX,sigX2,varargin)
%C_REG_EXP Online update coeffs. for regression with positive exponential noise.
%
%	Description
%
%	[LOGLIK,K1,K2] = C_REG_EXP(LIKPAR,Y,MX,VARX) - returns the
%	coefficients for the online update of the GP regression with single-
%	sided noise.
%
%	Parameters:
%
%	 LIKPAR  - the noise parameter (LAMBDA=NET.LIKPAR).
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
%	See also
%	OGP, OGPTRAIN, OGPFWD, ERR_ABS
%

%	Copyright (c) Lehel Csato (2001-2004)

varargout = {[]};

% checking the validity of the variance for the GP marginal
tolK2 = 1e-8;
if sigX2 < tolK2;
  warning(['Pred. distr. var.: ' num2str(sigX2)]);
  sigX2(find(sigX2<tolK2)) = tolK2;
end;

sqSig = sqrt(sigX2);
s     = sqSig./likpar;   % \lambda*sigma
y     = (mX-y)./sqSig;

% where the approximations should start
thresh = 6;
iAs    = find(y+s>thresh);
iN     = setdiff([1:length(y)],iAs);

tt            = 1./(y(iAs)+s(iAs)).^2;
l_f           = (1 - (5/2 - (37./3 - 353/4*tt).*tt).*tt).*tt;
logLik(iAs,1) = -(y(iAs).^2+log(2*pi))./2 - log(y(iAs)+s(iAs)) + l_f;

logLik(iN,1)  = s(iN).*(s(iN)/2+y(iN)) + ...
    log((1-erf((y(iN)+s(iN))./sqrt(2)))./2);

logLik    = logLik - log(likpar);

if nargout==1; return; end;   % STOP if nothing else is wanted

% go on with calculation
dt(iAs,1) = - (y(iAs)+s(iAs)).*exp(l_f);
dt(iN,1)  = - exp(-(y(iN)+s(iN)).^2/2)./ ...
    (1-erf((y(iN)+s(iN))/sqrt(2))).*2./sqrt(2*pi);

dt        = dt./sqSig;

K1        = 1./likpar + dt;

K2        = - dt.*( (y+s)./sqSig + dt);

% adjusting K2 - to make the INFERENCE alg. more stable
tolL  = 100;
indBad = find(-K2./(1 + K2.*sigX2)>tolL);
if length(indBad);
  K2(indBad) = - 1./(sigX2(indBad) + 1./tolL);
end;

indBad = find(-K2<1./tolL.^2);
if length(indBad);
  K2(indBad) = - 1./tolL.^2;
end;
