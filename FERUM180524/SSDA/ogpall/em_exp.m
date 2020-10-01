function [newLikPar] = em_exp(oldLikPar,y,cavM,cavV,postM,postV)
%EM_EXP	Recomputes the likelihood parameter for pos.-exp. noise.
%
%	Description
%
%	[NEWLIK] = EM_EXP(OLDLIKPAR,Y,CAVM,CAVV,POSTM,POSTV) -the M-step in
%	lik.par. estimation assuming one-sided exponential noise. The data is
%	assumed to be independent, and the EM approximation to the new value
%	is
%
%	newlik  = sum(avg(Theta(y-f),postM,postV))/N;
%
%	where THETA is a function that sets the negative values of its
%	argument to zero and AVG is the average operator with respect to a
%	Gaussian.
%
%	Parameters:
%
%	 OLDLIKPAR  - the preceeding value of the lik.-par..
%
%	 Y  - desired output.
%
%	 CAVM,CAVV  - (cavity) mean and variance of the GP.
%
%	 POSTM,POSTV  - mean and variance of the posterior process.
%	NEWLIK  - the new parameter.
%
%	See also
%	DEMOGP_REG, OGP, OGPTRAIN, C_REG_EXP
%

%	Copyright (c) Lehel Csato (2001-2004)

% checking the validity of the variance for the GP marginal
tolK2  = 1e-8;
bCav = find(cavV < tolK2);
bPost= find(postV < tolK2);
if length([bCav bPost]);
  warning('Variance TOO SMALL!');
  cavV(bCav) = tolK2;
  postV(bPost)=tolK2;
end;

% computing auxiliary quantities
s         = sqrt(postV);
x         = - (postM-y)./s/sqrt(2);
l         = ( s.*x.*(1+erf(x)) + exp(-x.^2)./sqrt(pi))./sqrt(2);
newLikPar = sum(l,1)./size(postV,1);

% implementing the momentum
lP = 0.5;
newLikPar = lP*newLikPar + (1-lP)*oldLikPar;
