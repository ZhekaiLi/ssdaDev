function [newLikPar] = em_lapl(oldLikPar,y,cavM,cavV,postM,postV)
%EM_LAPL Recomputes the likelihood parameter for exponential noise.
%
%	Description
%
%	[NEWLIK] = EM_LAPL(OLDLIKPAR,Y,CAVM,CAVV,POSTM,POSTV) - the M-step in
%	lik.par. estimation of the coefficeient of the Laplace noise. The
%	data is assumed to be independent, thus an analytic expression is
%	found (see MATLAB code).
%
%	Parameters:
%
%	 OLDLIKPAR  - preceeding value of the likelihood parameter.
%
%	 Y  - training output.
%
%	 CAVM,CAVV  - (cavity) mean and variance of the GP.
%
%	 POSTM,POSTV  - mean and variance of the posterior process.
%
%	 NEWLIK  - the new parameter.
%
%	See also
%	DEMOGP_REG, OGP, OGPTRAIN, C_REG_LAPL
%

%	Copyright (c) Lehel Csato (2001-2004)

% checking the validity of the variance for the GP marginal
tolK2  = 1e-8;
indBad = find(cavV < tolK2);
if length(indBad);
  warning(['Pred. distr. var.: ' num2str(min(cavV))]);
  cavV(indBad) = tolK2;
end;

s         = sqrt(abs(postV)); % making sure cannot be negative
x         = (postM-y)./s/sqrt(2);
l         = sqrt(2).*(s.*x.*erf(x) + exp(-x.^2)./sqrt(pi)) + 10*eps;
newLikPar = sum(l,1)./size(postV,1);

% implementing the momentum - THOUGHT to keep algorithm more stable
lP = 0.5;
newLikPar = lP*newLikPar + (1-lP)*oldLikPar;

if newLikPar<1/100;
  newLikPar = 1/100;
elseif newLikPar>100;
  newLikpar = 100;
end;
