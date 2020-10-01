function [allLog] = logp_l(likpar,y,mX,sigX2)
%LOGP_L	Computes the log-predictive probability for Gaussian noise
%
%	Description
%
%	[ALLLOG] = LOGP_L(LIKPAR,Y,MX,SIGX2) - returns the sum of log-
%	predictive probabilities of the inputs -- assuming additive Gaussian
%	noise.
%
%	It assumes that the test samples are independent, thus ignoring any
%	possible non-diagonal part in the joint posterior distribution.
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

% the average log-likelihood of the data for Laplace noise
logP = c_reg_lapl(likpar,y,mX,sigX2);
allLog = sum(logP)./length(y);
