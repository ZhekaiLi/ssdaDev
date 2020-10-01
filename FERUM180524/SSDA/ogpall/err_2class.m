function [errAvg, errInd] = err_2class(likpar,y,mX,sigX2)
%ERR_2CLASS Computes the binary classification error.
%
%	Description
%
%	[ERRAVG] = ERR_2CLASS(LIKPAR,Y,MX,SIGX2) - returns the percentage of
%	miclassified items.
%
%	Parameters:
%
%	 LIKPAR  - the likelihood parameter (not used here).
%
%	 Y  - training outputs.
%
%	 MX  - mean of the GP marginal at the new input.
%
%	 ALLERR  - percentage of misclassified inputs
%
%	[ERRAVG, ERRIND] = ERR_2CLASS(LIKPAR,Y,MX,SIGX2) - returns the
%	indicators for individual training data.
%
%	See also
%	DEMOGP_CLASS, DEMOGP_GUI, OGP, C_CLASS_BIN
%

%	Copyright (c) Lehel Csato (2001-2004)

% Checking whether the classification is binary
if ~(size(mX,2) == 1 & find(abs(y)==1)==length(y));
  error(['Function ''err_2class'' measures binary classification' ...
	 'error:\n GP output has to be one-dimensional and \n' ...
	 'desired values should be +/-1']);
end;

errInd = (y.*mX<0);
errAvg = mean(errInd)*100;
