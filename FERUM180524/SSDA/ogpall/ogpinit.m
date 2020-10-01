function ogpinit(likaddr, likpar,likoptfn);
%OGPINIT Likelihood initialisation for the Online Gaussian Process structure.
%
%	Description
%
%	OGPINIT(LIKFN, LIKPAR) uses the GLOBAL structure NET of an online
%	Gaussian Process data structure NET and sets the function (NET.LIKFN)
%	that returns the online update coefficients.  The vector NET.LIKPAR
%	sets the additional parameters to this function.
%
%	With three arguments: OGPINIT(LIKFN,LIKPAR,LIKOPTFN), the third
%	argument is the address of the optimiser to the likelihood parameter.
%	The optimiser assumes that there exists an approximation to the
%	posterior process and the values of the posterior mean and variance
%	at the training locations are available.  The optimisation of the
%	likelihood parameters is independent of the optimisation of the
%	covariance kernel parameters.
%
%	If the third argument in calling OGPINIT or the field NET.LIKOPTFN is
%	empty then there is no likelihood parameter optimisation in the
%	function OGPTRAIN.
%
%	Below there is a description of the functions and parameters involved
%	in likelihood optimsation.
%
%	The address in LIKFN is to a function that computes the coefficients
%	for the online learning.
%	[loglik, q, r] = likfn(likpar, y, mu, var, mu_p, var_p);
%	 with the following parameters:   LIKPAR  - the likelihood parameter
%	(see eg. C_REG_GAUSS).
%
%	 Y  - the training (noisy) outputs at the input.
%
%	 MU, VAR  - the predicted mean and variance (i.e. the   statistics of
%	the marginal GP at the input).     MU_P, VAR_P  - the values of the
%	prior mean and variance   before the addition of the new data -
%	called cavity means and variances.   Providing these values is
%	optional and it is returned by OGPPOST within   the structure GPOPT.
%
%	 LOGLIK  - the value of the log-average.
%
%	 Q, R  - the update coeficients to the online learning   algorithm -
%	vector and matrix of size NOUT (see C_REG_GAUSS).
%
%
%	The function LIKOPTFN has the following structure:
%
%	newlikpar = likoptfn(oldlikpar, y, cavM, cavV, postM, postV);
%	 where the input and output parameters are the following:
%
%	 NEWLIKPAR, OLDLIKPAR  - the new and old values of the   likelihood
%	parameters.
%
%	 Y  - the vector of training inputs
%
%	 CAVM, CAVV  - vector of prior means and variances   corresponding to
%	the training locations !! AND !! with the contribution of   the
%	current input removed (cavity parameters).
%
%	 POSTM, POSTV  - (optional) vector of posteior means and   variances
%	which is sometimes needed, e.g. if one wants to use the EM
%	algorithm.
%
%
%	See also
%	OGP, OGPTRAIN, C_REG_GAUSS, C_REG_EXP, EM_GAUSS, EM_LAPL
%

%	Copyright (c) Lehel Csato (2001-2004)

global net;

% consistency checking
if nargin<2;
  error('Wrong parameter to function OGPINIT');
end;

net.likaddr  = likaddr;

if nargin>=2;		      % initial value for LIK.PAR is given
  net.likpar = likpar;
else
  net.likpar = 0;
end;

if nargin>=3;		      % function to optimise LIK.PAR
  net.likoptfn = likoptfn;
else
  net.likoptfn = [];
end;
