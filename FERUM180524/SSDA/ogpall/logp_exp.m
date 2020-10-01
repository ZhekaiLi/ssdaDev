function [allLog] = logp_exp(likpar,y,mX,sigX2)
%LOGP_EXP Computes the log-predictive probability for positive exponential noise
%
%	Description
%
%	[ALLLOG] = LOGP_EXP(LIKPAR,Y,MX,SIGX2) - returns the sum of log-
%	predictive probabilities of the inputs -- assuming additive positive
%	exponential noise.  The samples are treated independently.
%
%	Parameters:
%
%	 LIKPAR  - the likelihood parameters.
%
%	 Y  - desired output.
%
%	 MX  - mean of the GP marginal at the new input.
%
%	 ALLLOG  - the average error.
%
%	See also
%	OGP, OGPTRAIN, DEMOGP_REG, ERR_ABS
%

%	Copyright (c) Lehel Csato (2001-2004)

% the average log-likelihood of the data - POSEXP NOISE MODEL!!!
[logP] = c_reg_exp(likpar,y,mX,sigX2);
allLog = sum(logP)./length(y);
