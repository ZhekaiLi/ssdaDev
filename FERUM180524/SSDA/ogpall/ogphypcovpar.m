function ogphypcovpar(hyplambda, hypmean)
%OGPHYPCOVPAR Initialises the Sparse Gaussian Process hyperparameters.
%
%	Description
%
%	OGPHYPCOVPAR(HYPLAMBDA) initialises the hyperparameter vector of the
%	Sparse Gaussian Process structure NET, where HYPLAMBDA specifies the
%	precision (inverse variance) of the parameters - where the parameter
%	order is given by the function OGPPAK.  Since several parameters have
%	to be positive, OGPPAK usually contains the log-transformed
%	parameters, forcing positivity.
%
%	If a single parameter is given, zero mean is assumed to the
%	(log)parameters except for the input scales which are initialised to
%	-4.605-LOG(2*NIN) which supposes a preference for lengthscales of
%	order 100 increasing with the number of dimensions.
%
%	If unnormalised inputs are used, then it is advised to have the mean
%	values initialised to values proportional to negative log-variance.
%	The inference on the hyper-parameter level is ML2 inference, thus
%	vague priors are beneficial for the stability of the algorithm.
%
%	A manual setup of the hyperparameters is recommended only if some
%	specific knowledge about the data is known.
%
%	The function sets the fields of the GLOBAL variable NET.
%
%	See also
%	OGP, OGPPAK, OGPCOVARP, OGPTRAIN
%

%	Copyright (c) Lehel Csato (2001-2004)

global net;

% checking parameters
if nargin==0;
  hyplambda = 0;
end;

hpL = size(ogppak,2);
if size(hyplambda,1)>1;
  error('Single line allowed for hyperparameters');
elseif length(hyplambda)==1 & hyplambda;
  hyplambda = hyplambda*ones(1,hpL);
end;
net.hyplambda = hyplambda;

if nargin==1;
  hypmean = zeros(1,hpL);
  hypmean(1) = -1000;	      % ZERO BIAS
  % Initialising the prior means for log-lengthscales
  hypmean(2:(net.nin+1)) = - min(1.5 + log(2*net.nin) , 10);
  % prior for SIGNAL amplitude
  hypmean(2+net.nin)     = log(1);
  % prior for signal order
  hypmean(3+net.nin)     = log(2);
end;
net.hypmean = hypmean;
