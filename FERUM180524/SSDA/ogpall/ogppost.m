function ogppost(xTrain,yTrain)
%OGPPOST Calculation of the sparse posterior
%
%	Description
%
%	OGPPOST(X_TRAIN,Y_TRAIN,ITS) - considers the GLOBAL (prior) Gaussian
%	Process structure NET together with the training data
%	[X_TRAIN,Y_TRAIN] and sets the (posterior) Gaussian Process structure
%	NET using the training data within the online learning algorithm.
%	The number of sweeps through the dataset is ITS, if unspecified, its
%	default value is 1.
%
%	In Sparse online learning, by default the BV set changes, which
%	implies significant matrix manipulation.  If we want the BV set to be
%	constant, we need to set the indicator ISBVFIXED (e.g.
%	NET.ISBVFIXED=1).  This makes the code faster, but it also implies
%	that one has to set the BV set ``manually'' (using OGPEMPTYBV).
%	Alternatively one should first make an iteration using the default BV
%	set selection mechanism and then set NET.ISBVFIXED=1.
%
%	The functions OGPOST and OGPTRAIN use other two GLOBAL structures
%	OGPTRAIN and EP.  Setting the fields of these structures influences
%	the way the posterior is computed.
%
%	OGPPOST(X_TRAIN,Y_TRAIN) - uses the GLOBAL structure GPOPT to
%	transmit/receive additional variables from/to the training, like eg.
%	the TAP/EP extension of the online learning, the test- or training
%	errors and the various options to perform hyperparameter optimisation
%	(see DEFOPTIONS for details).
%
%	Two sets of hyperparameters are distinguished: those related to the
%	kernel function and those characterising the likelihood. The
%	initialisation of the kernel parameters is done in OGP and those of
%	the likelihood function in OGPINIT (see the respective functions).
%	Optimisation of the hyperparameters is implemented within the
%	function OGPTRAIN.
%
%	See also
%	OGP, OGPINIT, OGPTRAIN, OGPSTEP_FULL, OGPSTEP_SP, DEFOPTIONS
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

%%%%%%%%%%%%%%%%%%%%% ARGUMENT CHECKING %%%%%%%%%%%%%%%%%%%%%%%%%
nTrain    = size(xTrain,1);

if isfield(gpopt.postopt,'isep') && gpopt.postopt.isep;
  if isempty(ep) || ~isfield(ep,'lamP') || isempty(ep.lamP);
    % initialise the EP parameters
    ep.lamP  = sparse(nTrain*net.nout,nTrain*net.nout);
    ep.aP    = zeros(nTrain*net.nout,1);
    ep.projP = zeros(nTrain*net.nout,size(net.BV,1)*net.nout);
    ep.logZ  = zeros(nTrain,1);
    ep.X     = xTrain;
  end;
end;

%%%%%%%%%%%%%%%%%%%%% THE ONLINE ITERATIONS %%%%%%%%%%%%%%%%%%%%%
allSteps = nTrain*gpopt.postopt.itn;
for iSweep = 1:gpopt.postopt.itn;
  if gpopt.postopt.shuffle;   %	inputs in a random order
    tPerm = randperm(nTrain);
  else
    tPerm = 1:nTrain;
  end;
  for pInd=1:nTrain;	      % considers the data items
    % a 'long index'
    cInd = tPerm(pInd);
    lInd = pInd + nTrain * (iSweep-1);
    % get current data
    iKX   = ogpcovarp(net.BV,xTrain(cInd,:));
    ikk   = ogpcovarp(xTrain(cInd,:));
    if ~strcmp(upper(net.outtype),'NONE');
      iCY   = yTrain(cInd,:);
    else
      iCY   = [];
    end;
    %%%%%%%%%%% REMOVING THE CURRENT SITE %%%%%%%%%%%%%
    nBV  = size(net.BV,1);
    indS = [1:net.nout]  + (cInd-1)*net.nout;
%     if gpopt.postopt.isrm  % flag to remove BV or not
    if gpopt.postopt.isep;
      if abs(det(ep.lamP(indS,indS)))>1e-10;
	teK = ep.projP(indS,:) * net.KB;
	tMe = (eye(nBV*net.nout) + net.C*net.KB) * ...
	      ep.projP(indS,:)';
	tV  = ep.lamP(indS,indS);
	tV  = (eye(net.nout) - tV*teK*tMe )\tV;
	tV  = (tV + tV')./2;
	net.w = net.w + tMe * tV*(teK*net.w - ep.aP(indS));
	net.C = net.C + tMe * tV * tMe';
      end;
    end;
%     end
    %%%%%%%%%%% CAVITY MEAN/VARIANCE %%%%%%%%%%%%%%%%%%
    cM  = iKX'*net.w;
    cV  = ikk + iKX'*net.C*iKX;
    cV  = (cV + cV')./2;

    %%%%%%%%%%% CHECKING VARIANCE %%%%%%%%%%%%%%%%%%%%%
    if net.nout>1;
      [eVec,eVal] = eig(cV);
      eVal        = diag(eVal);
      iWrong      = find(eVal<1e-12);
      if length(iWrong);
	fprintf('Correcting prior!\n')
	eVal(iWrong) = 1e-12;
      end;
      eVal        = diag(eVal);
      cV          = eVec*eVal*eVec';
    else 
      cV = max(cV,1e-12);
    end;

    %%%%%%%%%%% COMPUTING THE PRIOR MEAN %%%%%%%%%%%%%%
    if ~isnumeric(net.prmean);
      pM = feval(net.prmean,xTrain(cInd,:),net.prmeanp);
    else
      pM = zeros(net.nout,1);
    end;
    %%%%%%%%%%% STORING - FOR POST-PROCESSING %%%%%%%%%
    ep.cavM(indS,:) = cM+pM;
    ep.cavV(indS,:) = cV;

    %%%%%%%%%%% COMPUTING THE ONLINE COEFFS %%%%%%%%%%%
    [logEvid, K1, K2] = ...
	feval(net.likaddr,net.likpar,iCY',cM+pM,cV,pM,ikk);

    %%%%%%%%%%% GEOMETRIC COMPONENTS %%%%%%%%%%%%%%%%%%
    hatE  = net.KB\iKX;
    gamma = ikk - iKX'*hatE;
    %%%%%%%%%%% ADJUSTING TO GIVE VALID VERS %%%%%%%%%%
    % the last parameter specifies how negative the VAR. can be
    [K1, K2]  = ogpparadj(K1,K2,cM+pM,cV,5e3*ikk,1e-10);

    %%%%%%%%%%% OGP PARAMETER UPDATE %%%%%%%%%%%%%%%%%%
    if (det(gamma) < net.thresh*det(ikk)) || net.isBVfixed;
      ogpstep_sp(cInd, iKX, K1, K2, gamma, hatE);
      sp_step = 1;
    else
      ogpstep_full(xTrain(cInd,:), iKX, ikk, K1, K2, ...
			 gamma, hatE);
      sp_step = 0;
    end;
    if gpopt.postopt.isep;    % updating EP parameters
      ogpstep_ep(sp_step,cInd,K1,K2, ...
		 gamma,hatE,cM,cV,logEvid);
    end;

    %%%%%%%%%%% REMOVING EXTRA BV.S %%%%%%%%%%%%%%%%%%%
    if ~net.isBVfixed;
      while 1;
	nBV = size(net.BV,1);
	[iMin, allScore] = ogpbvmin;

	if (nBV <= net.maxBV); break; end;
	if gpopt.postopt.isep;
	  ogpdelbv(iMin);
	else
	  ogpdelbv(iMin);
	end;
      end;
    end;
    %%%%%%%%%%% REMOVAL BASED ON GEOMETRY %%%%%%%%%%%%%
    while nBV;
      [iMin,gScore] = ogpbvmin(1); % only geometry

      if gScore(iMin) >= net.thresh/1000; break; end;

      ogpdelbv(iMin);
      nBV = size(net.BV,1);
      fprintf('.');
    end;

    %%%%%%%%%%% PROCESSING GPOPT %%%%%%%%%%%%%%%%%%%%%%
    if gpopt.pavg | gpopt.ptest | gpopt.ptrain;
      if gpopt.pavg;	      % log-avg for a single data
	gpopt.logavg(lInd) = logEvid;
      end;
      % errors
      if ~mod(allSteps-lInd,gpopt.freq);
	tInd = 1+ max(size(gpopt.testerror,2),size(gpopt.trainerror,2));
	% test error
	if gpopt.ptest;
	  [meanT varT] = ogpfwd(gpopt.xtest);
	  for iErr = 1:length(gpopt.erraddr);
	    gpopt.testerror(iErr,tInd) = feval( ...
		gpopt.erraddr{iErr},net.likpar,gpopt.ytest,meanT,varT);
	  end;
	end;
	% training error
	if gpopt.ptrain;
	  [meanT varT] = ogpfwd(xTrain);
	  for iErr = 1:length(gpopt.erraddr);
	    gpopt.trainerror(iErr,tInd) = feval( ...
		gpopt.erraddr{iErr},net.likpar,yTrain,meanT,varT);
	  end;
	end;
	% displaying results if needed
	if gpopt.disperr & ...
	      (gpopt.ptrain | gpopt.ptest | gpopt.pavg);
	  fprintf('Step %4d, %4d BVs -',lInd,size(net.BV,1));
	  if gpopt.ptrain;
	    fprintf(' trnE:');
	    for iErr = 1:length(gpopt.erraddr);
	      fprintf(' %1d:%5.3f;',iErr,gpopt.trainerror(iErr,tInd));
	    end;
	  end;
	  if gpopt.ptest;
	    fprintf(' tstE: ');
	    for iErr = 1:length(gpopt.erraddr);
	      fprintf('%1d: %5.3f;',iErr,gpopt.testerror(iErr,tInd));
	    end;
	  end;
	  if gpopt.pavg;
	    fprintf('  log-avg.: %5.3f;',logEvid);
	  end;
	  fprintf('\n');
	end;
      end;
    end;		      % END processing GPOPT
  end;			      % END inner cycle
end;			      % END outer cycle
