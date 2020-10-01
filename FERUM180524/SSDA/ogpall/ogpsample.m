function [samp] = ogpsample(x, approx, kxx, kxbv,rseed)
%OGPSAMPLE Generates samples from a Gaussian process.
%
%	Description
%
%	OGPSAMPLE(X) takes a GLOBAL Spase OGP data structure NET together
%	with a matrix X of input vectors and returns a matrix of realisations
%	of the marginal Gaussian distribution at the inputs X
%
%	SAMP = OGPSAMPLE(NET, X, APPROX),with the binary indicator APPROX set
%	to nonzero,samples from a restricted subspace by inverting only the
%	covariance matrix correpsonding to the BV set and ignoring the Schuur
%	complement of KXX with respect to KBV.
%
%	SAMP = OGPSAMPLE(NET, X_INPUTS, APPROX, KXX, KXBV) uses the pre-
%	computed prior covariance matrices KXX and KXBV to generate the
%	samples. This increases efficiency if several calls to OGPSAMPLE are
%	made.
%
%	The time required for this sampling procedure is cubic with respect
%	to the size of the inputsdue to the inversion of the covariance
%	matrix.  Use it with caution.
%
%	If APPROX is used, then there are ``only'' as many random variables
%	generated as many BV points are in the Gaussian process structure
%	NET, meaning that matrix inversion of the size of the BV set is to be
%	performed.
%
%	SAMP = OGPSAMPLE(NET, X, APPROX, KXX, KXBV, RSEED) includes the
%	option of specifying the seed for the random number generator. Used
%	to generate the same random numbers in repeated samples -- for
%	different kernels.  If the matrices are not computed, then empty
%	matrices [] should be used instead.
%
%	See also
%	OGP, OGPINIT, OGPFWD
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% consistency checking
errstring = consist(net, 'ogp', x);
if ~isempty(errstring);
  error(errstring);
end

bvDim = size(net.BV,1) * net.nout;
xDim  = size(x,1) * net.nout;
% checking extra arguments
if nargin<3;
  approx = 0;		      % default value for approx
end;
if nargin<4		      % computing covariance
  kxbv = ogpcovarp(x, net.BV);
  kxx  = ogpcovarp(x, x);
  kxx  = (kxx + kxx')./2;
else
  if any(size(kxx)~=[xDim xDim]);
    kxx  = ogpcovarp(x, x);
    kxx  = (kxx + kxx')./2;
  end;
  if any(size(kxbv)~=[xDim bvDim]);
    kxbv = ogpcovarp(x, net.BV);
  end;
end;
% adding the random seed choice
if nargin==6;		      % SEED is the last argument
  randn('state',rseed);
end;

% COVARIANCE
if ~approx;
  Cov  = kxx + kxbv*net.C*kxbv';
  Cov  = (Cov + Cov')./2;
  [vCov, dCov] = eig(Cov);
  samp = randn(xDim,1);
else
  Cov  = net.KBinv + net.C;
  Cov  = (Cov + Cov')./2;
  [vCov, dCov] = eig(Cov);
  vCov = kxbv*vCov;
  samp = randn(bvDim,1);
end;
  dCov = sqrt(abs(diag(dCov)));
  dCov = diag(dCov);
% GENERATING data
samp = kxbv*net.w + vCov*dCov * samp;
% RESHAPING data
samp = reshape(samp,[size(x,1) net.nout]);
