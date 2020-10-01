function [fopt]=defoptions(varargin);
%DEFOPTIONS sets the ``gpopt'' structure to default values.
%
%	Description
%
%	GPOPT = DEFOPTIONS initialises GPOPT, the structure that controls the
%	training procedure OGPTRAIN. Below is a description of the GPOPT
%	structure.
%
%	The fields in GPOPT are divided into options related to the
%	approximation of the posterior process (GPOPT.POSTOPT); options
%	related to the optimisation of the covariance function parameters
%	(GPOPT.COVOPT); and parameters returned by the optimisation procedure
%	like test or training errors during training or the value of the
%	predictive log-likelihood. The structure GPOPT also has options to
%	compute the the test/training error, the marginal likelihood for the
%	test/training data, etc.
%
%	Fields of GPOPT:
%
%	 POSTOPT  - substructure with parameters related to the computation
%	of the posterior process - keeping the covariance and likelihood
%	parameters fixed.
%
%	 COVOPT  - substructure grouping parameters for the optimisation of
%	the covariance parameters.
%
%	 PAVG  - boolean indicator for storing (or NOT) the log-averages. If
%	nonzero, LOGAVG stores the sequence of log-averages for each training
%	input.
%
%	 DISPERR  - boolean indicator whether to display (NOT) the errors
%	during training.
%
%	 ERRADDR  - address of function to compute the test error.  It should
%	have four inputs: NET, the desired outputs Y, the predictive means
%	M_X, and the predictive variances VAR_X.  The function ERRADDR
%	returns a (user-specified) measure of error. If there is no function
%	given, the weighted quadratic error (implemented in ERR_MSE) is used.
%	See also this function on how to implement new error functions.
%
%	 PTEST  - indicator to store test errors.  Evaluating test error can
%	be expensive, GPOPT.FREQ specifies the delays between successive test
%	error computation (0=1, i.e. test error for each online step).
%
%	 X_TEST,Y_TEST  - the test inputs and outputs.
%
%	 TESTERROR  -  the returned test errors.
%
%	 PTRAIN  - indicator whether to compute or not the training errors.
%	If this value is nonzero then, similarly to computing test errors,
%	the training errors are computed GPOPT.FREQ-th step.
%
%	 TRAINERROR  - the returned training errors.
%
%	The structure GPOPT.POSTOPT stores options driving the computation of
%	the posterior process:
%
%	 POSTOPT.ITN  - number of online sweeps through the data (default 1).
%
%	 POSTOPT.SHUFFLE  - if nonzero (by default), then the inputs are
%	shuffled at each iteration, this is an attempt to make the posterior
%	independent of the data ordering.
%
%	 POSTOPT.ISEP  - indicator whether to perform the TAP/EP learning
%	procedure or not.  This requires additional values to be kept for
%	further processing.
%
%	 POSTOPT.FIXITN  - keeps the basis vectors fixed and performs the
%	TAP/EP iteration.  Thus one source of fluctuations is eliminated, and
%	the TAP/EP parameters become stable.
%
%	If GPOPT.POSTOPT.ISEP is set to nonzero, then the inference uses the
%	TAP/EP iterative approach, which is more time-consuming and also
%	requires additional additional information to be stored.  These
%	values are stored in the substructure GPOPT.EP using the following
%	fields:
%
%	 EP.X  - location of training inputs.
%
%	 EP.AP  - mean of the site distribution corresponding to the
%	likelihood.
%
%	 EP.LAMP  - site variances corresponding to the likelihood.
%
%	 EP.PROJP  - coefficients of the projection.
%
%	The substructure GPOPT.COVOPT contains the fields related to the
%	optimisation of the covariance function parameters.  The optimisation
%	relies on the NETLAB optimisation routines, the default is
%	'CONJGRAD', this is the string stored in FNOPT. Additional options to
%	the respective optimisation routine are provided via the field OPT.
%	If this is a scalar, then is specifies the number of iterations,
%	otherwise it has to conform to the NETLAB specifications (see NETOPT
%	from the NETLAB package).
%
%	See also
%	OGP, OGPTRAIN, EM_GAUSS
%

%	Copyright (c) Lehel Csato (2001-2004)

% default values
fopt = [];
fopt.postopt         = [];    % options for computing the posterior GP
fopt.postopt.itn     = 1;     % (def): single iter.
fopt.postopt.fixitn  = 0;     % (def): NO extra iteration with fixed BV
fopt.postopt.shuffle = 1;     % (def): shuffle data set
fopt.postopt.isep    = 0;     % (def): no EP-step

fopt.disperr = 0;     	      % (def): no DISPLAY
fopt.pavg    = 0;     	      % (def): no output
fopt.ptest   = 0;     	      % (def): no test set
fopt.ptrain  = 0;     	      % (def): no train error
fopt.testerror = [];
fopt.trainerror= [];
fopt.erraddr = {};    	      % NO ERROR COMPUTATION

fopt.covopt = [];             % options for cov. par. optimisation
fopt.covopt.fnopt = 'scg';% function used to find the minimum
fopt.covopt.opt = [];	      % number of iterations - default - 2x(par#)
