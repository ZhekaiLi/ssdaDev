function [allLog] = logp_g(likpar,y,mX,sigX2)
%LOGP_G	Computes the log-predictive probability for Gaussian likelihood
%
%	Description
%
%	[ALLLOG] = LOGP_G(LIKPAR,Y,MX,SIGX2) - returns the average log-
%	predictive probabilities of the inputs -- assuming additive Gaussian
%	noise.
%
%	It assumes that the test samples are independent, thus ignoring any
%	possible non-diagonal part in the joint posterior distribution.
%
%	Parameters:
%
%	 LIKPAR  - the current likelihood parameters.
%
%	 Y  - desired output.
%
%	 MX  - mean of the GP marginal at the new input.
%
%	 ALLLOG  - the average error.
%
%	See also
%	OGP, OGPTRAIN, DEMOGP_REG, ERR_MSE
%

%	Copyright (c) Lehel Csato (2001-2004)

% the average log-likelihood of the data - GAUSSIAN noise
tVal   = sigX2 + likpar(1);  
allLog = -log(2*pi) - sum(log(tVal));
tVal   = (y-mX).^2./tVal;
allLog = (allLog - sum(tVal))./2./length(y);
