function [newLikPar] = em_gauss(oldLikPar,y,cavM,cavV,postM,postV)
%EM_GAUSS Recomputes the (Gaussian) noise variance. 
%
%	Description
%
%	[NEWLIKPAR] = EM_GAUSS(OLDLIKPAR,Y,CAVM,CAVV,POSTM,POSTV) -
%	implements the M-step in the GP inference assuming Gaussian additive
%	noise. The data is assumed to be factorising, thus we have a direct
%	result as
%	newLikPar  = sum(postV + (y-postM).^2) / N
%	 where N is the number of training samples.
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
%	NEWLIKPAR  - the new noise variance.
%
%	See also
%	OGPTRAIN, DEMOGP_REG, C_REG_GAUSS, G_L_GAUSS
%

%	Copyright (c) Lehel Csato (2001-2004)

% checking the validity of the variance for the GP marginal
tolK2  = 1e-8;
bCav = find(cavV < tolK2);
bPost= find(postV < tolK2);
if length([bCav; bPost]);
  warning('Variance TOO SMALL!');
  cavV(bCav) = tolK2;
  postV(bPost)=tolK2;
end;

%% THE EM-procedure which relies on an approximation to the marginal lik.
newLikPar = sum(postV + (y-postM).^2,1)/size(postV,1);

% momentum parameter needed since the estimation of LIKPAR is more often
% -- hope to prevent overfitting.
lP = 0.5;
newLikPar = lP*newLikPar + (1-lP)*oldLikPar;
