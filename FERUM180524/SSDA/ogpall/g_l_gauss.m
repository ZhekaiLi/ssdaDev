function [newLikPar] = g_l_gauss(oldLikPar,y,cavM,cavV,postM,postV)
%G_L_GAUSS Recomputes the variance of the Gaussian noise using gradients
%
%	Description
%
%	[NEWLIKPAR] = G_L_GAUSS(OLDLIKPAR,Y,CAVM,CAVV,POSTM,POSTV) -
%	implements a conjugate gradient (SCG) algorithm for adjusting the
%	noise variance in the Gaussian noise model.  The gradient steps are
%	based on the approximation to the marginal likelihood (evidence)
%	using the cavity distributions obtained with the TAP/EP framework.
%
%	This version uses the cavity means and variances of the GP
%	approximation - unlike the EM-based method (EM_GAUSS) which uses the
%	posterior GP.
%
%	The algorithm employs the scaled conjugate gradient (SCG) method from
%	the NETLAB toolbox.
%
%	A momentum term is introduced whih is aimed to stabilise the
%	algorithm.
%
%	Parameters:
%
%	 OLDLIKPAR  - the old (previous) value of the Gaussian noise
%	parameter.
%
%	 Y  - desired output.
%
%	 CAVM,CAVV  - (cavity) mean and variance of the GP.
%
%	 POSTM,POSTV  - mean and variance of the posterior process.
%	NEWLIKPAR  - the new noise variance.
%
%	See also
%	OGPTRAIN, DEMOGP_REG, C_REG_GAUSS
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

gradLik = inline('x*sum((1 - mm2./(x.^2 +sigX2))./(x.^2 + sigX2))', ...
		 'x', 'sigX2', 'mm2');
origLik = inline('sum((log(x.^2 + sigX2) + mm2./(x.^2 +sigX2)))./2', ...
		 'x', 'sigX2', 'mm2');

nnopt     = foptions;	      % preparing NETLAB optimisation
nnopt(1)  = -1;		      % ABSOLUTELY NO feedback!!
nnopt(14) = 10;		      % number of SCG steps
nnopt(9)  = 0;		      % FOR GRADCHECK!!!! NOT NEEDED!!!
stdNew    = scg(origLik,sqrt(oldLikPar),nnopt,gradLik,cavV,(y-cavM).^2);
newLikPar = stdNew^2;

% momentum parameter needed since the estimation of LIKPAR is more often
% -- hope to prevent overfitting.
lP = 0.7;
newLikPar = lP*newLikPar + (1-lP)*oldLikPar;
