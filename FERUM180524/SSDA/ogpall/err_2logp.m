function [logAvg, logInd] = err_2logp(likpar,y,mX,sigX2)
%ERR_2LOGP Computes the log-predictive probability of the labels.
%
%	Description
%
%	[LOGAVG] = ERR_2LOGP(LIKPAR,Y,MX,SIGX2) - average log-predictive
%	likelihood for the data-set.
%
%	Parameters:
%
%	 LIKPAR  - the likelihood parameters.
%
%	 Y  - desired output.
%
%	 MX  - mean of the GP marginal at the new input.
%
%	 ALLERR  - the mean-square error.
%
%	[LOGAVG, LOGIND] = ERR_2LOGP(LIKPAR,Y,MX,SIGX2) - additionally to the
%	average, it also returns the individual errors.
%
%	See also
%	DEMOGP_CLASS, DEMOGP_GUI, OGP, C_CLASS_BIN
%

%	Copyright (c) Lehel Csato (2001-2004)

% Checking whether the classification is binary
if ~(size(mX,2) == 1 & find(abs(y)==1)==length(y));
  error(['Function ''err_2logp'' measures binary classification' ...
	 'error:\nOGP output has to be one-dimensional and \n' ...
	 'test outputs should be +/-1']);
end;

mX     = y.*mX./sqrt(likpar(1) + sigX2);
logInd = zeros(size(mX));
% avoiding numerical errors: splitting the inputs...
iSm    = find(mX<-5);
iNn    = find(mX>=-5 & mX<6);
iLl    = find(mX>=6);

logInd(iNn) = - log((1 + erf(mX(iNn)./sqrt(2)))./2);
tt = mX(iSm).^(-2);
logInd(iSm) = .5*(log(2*pi) + mX(iSm).^2) + log(-mX(iSm)) - ...
    (1-(2.5-(37/3-353/4.*tt).*tt).*tt).*tt;
logInd(iLl) = 0;

logAvg = sum(logInd)./length(y);
