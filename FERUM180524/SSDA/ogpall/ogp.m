function ogp(nin, nout, covarfn, covarfnpar, prmean, prmeanp)
%OGP	Initialises the global net structure for the OGP toolbox.
%
%	Description
%
%	OGP(NIN, NOUT, COVARFN,COVPAR) creates a Gaussian Process model in
%	the GLOBAL structure NET which has NIN input and NOUT output
%	dimensions.
%
%	The string in the field COVARFN specifies the type of the covariance
%	function to be used.  The parameters to the covariance function are
%	given in COVPAR.  The available covariance functions are listed in
%	OGPCOVARF.
%
%	The function returns a data structure NET with the rest of parameters
%	set to zero. If COVPAR is not specified then the default values are
%	assigned to it (description of implemented kernel functions is in
%	OGPCOVARF).
%
%	OGP(NIN, NOUT, COVARFN, COVPAR, PRMEAN, PRMEANP) also sets a prior
%	mean function to the Gaussian process. The address in PRMEAN is the
%	function returning the prior means at locations CURRX. The field
%	PRMEANP contains optional parameters to PRMEAN.
%
%	The function PRMEAN has the structure
%
%	[meanVec] = prmean(x,prmeanp);
%
%	The likelihood function and its parameters can be set in OGPINIT,
%	similarly to PRMEAN (see OGPINIT).
%
%	Additional parameters which influence the computation of the
%	posterior GP and the re-calculation of the hyperparameters are
%	changed directly, using the designed fields of the structure NET.
%
%	The structure NET is global, thus there is no need for it to be
%	transmitted as a parameter.
%
%	See the online (HTML) documentation for details
%
%	See also
%	OGPINIT, OGPCOVARF, OGPPAK, OGPUNPAK, OGPFWD
%

%	Copyright (c) Lehel Csato (2001-2004)

global net;

net.type = 'ogp';
net.nin = nin;
net.nout = nout;

% Additional working parameters
net.thresh    = 1e-8;
net.maxBV     = 200;
net.KBinv     = zeros(0,0);
net.KB        = zeros(0,0);

% known covariance functions
covarfns = upper({'sqexp','ratquad','poly','linspline','fsin','sinc', ...
		  'polyexp', 'ou', 'user'});

if ~sum(strcmp(upper(covarfn), covarfns));
  error('Undefined activation function. Exiting.');
else
  net.covarfn = covarfn;
end

% some sort of default values would be welcome
net.kpar = covarfnpar;

% Initialise the GP/EP parameters
net.bias      = -1000;
net.hyplambda = 0;	      % prior inv.var. for Hyp.Pars
net.hypmean   = 0;	      % prior mean for Hyp.Pars

net.BV        = zeros(0,net.nin);
net.w         = zeros(0,1);
net.C         = zeros(0,0);
net.isBVfixed = 0;
net.likaddr   = [];	      % is initialised at opginit
net.likpar    = [];	      % is initialised at opginit
net.outtype   = 'DIRECT';     % type of the likelihood
if nargin>4;
  net.prmean  = prmean;	      % setting up the prior mean function
  net.prmeanp = prmeanp;      % and the afferent variables;
else
  net.prmean  = [];	      % NO prior mean function
end;

% DEALING with KL-projections the other way round: nonzero value
% indicates that the projection which keeps the canonical var.s (EP ones)
% constant, ZERO value indicates that the GP marginals at BV points are
% conserved!!! 
net.proj = 1;

% TEMPORARY variables
net.addEps = 1;
