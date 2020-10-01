function [absAvg, absInd] = err_abs(likpar,y,mX,sigX2)
%ERR_ABS Computes the absolute error.
%
%	Description
%
%	[ABSAVG] = ERR_ABS(LIKPAR,Y,MX,SIGX2) - returns the sum of absolute
%	errors.
%
%	The additional second output argument returns the individual absolute
%	errors:  [ABSAVG, ABSIND] = ERR_ABS(LIKPAR,Y,MX,SIGX2)
%
%	The input argument LIKPAR is not used.
%
%	Parameters:
%
%	 LIKPAR  - the Gaussian Process data structure.
%
%	 Y  - desired output.
%
%	 MX  - mean of the GP marginal at the new input.
%
%	 ABSAVG  - the average absolute error.
%
%	 ABSIND  - the vector of individual absolute errors.
%
%
%	See also
%	OGP, OGPTRAIN, DEMOGP_REG, C_REG_LAPL, C_REG_EXP
%

%	Copyright (c) Lehel Csato (2001-2004)

%% the average ABSOLUTE error
absInd = sum(abs(y-mX),2);
absAvg = mean(absInd);
