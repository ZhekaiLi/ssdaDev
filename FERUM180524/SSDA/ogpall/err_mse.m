function [mseAvg, mseInd] = err_mse(likpar,y,mX,sigX2)
%ERR_MSE Computes the mean-square error.
%
%	Description
%
%	[MSEAVG] = ERR_MSE(LIKPAR,Y,MX,SIGX2) - returns the mean-square error
%	between the training outputs Y and the posterior MX.
%
%	Called with two output arguments, the function returns the individual
%	squared distances as a second output argument.
%
%	The input argument LIKPAR is not used.
%
%	Parameters:
%
%	 NET  - the Gaussian Process data structure.
%
%	 Y  - desired output.
%
%	 MX  - mean of the GP marginal at the new input.
%
%	 ALLERR  - the mean-square error.
%
%	See also
%	DEMOGP_REG, OGP, C_REG_GAUSS
%

%	Copyright (c) Lehel Csato (2001-2004)

%% mean-square errors
mseInd = sum((y-mX).^2,2);
mseAvg = mean(mseInd);
