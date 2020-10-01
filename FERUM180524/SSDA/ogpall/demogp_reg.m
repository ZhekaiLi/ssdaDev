%DEMOGP_REG One-dimensional regression example with different noise models.
%
%	Description
%
%	DEMOGP_REG - script to approximate the noisy SINC function.
%
%	There is a possibility of changing the noise type from Gaussian
%	('GAUSS') to symmetric ('LAPLACE') and nonsymmetric ('POSEXP')
%	exponential distribution.
%
%	The script computes the test error and displays its average for
%	different BV set sizes.  The evolution of the average log-evidence is
%	also plotted.
%
%	This script can be used as copy-paste in own application.
%
%	See also
%	OGP, OGPTRAIN, C_REG_GAUSS, C_REG_EXP, C_REG_LAPL, SINCDATA
%

%	Copyright (c) Lehel Csato (2001-2004)

% clearing the workspace 
clear all;
global net gpopt ep;

% General information
allBV  = 400;	      % BV set sizes
totExp = 1;		      % number of experiments for EACH BV set
totHyp = 1;		  % Number of EM-steps in learning model hyperparameters
totSCG = 8;		      % number of SCG steps scaled conjugate gradient
nTest  = 201;
nTrain = 200;

% the type of and magnitude (or inverse magn.) of the additive noise.
nType  = 'gauss';	      % noise type: gauss, laplace, posexp
sig02  = 0.;		      % noise variance

% TEST/VISUALISATION data
[xTest,yTest] = sincdata(nType,nTest,0,[],1);

figure(4343);
indBV = 0;
for maxB = allBV;
  indBV = indBV + 1;
  fprintf('\n %4d :\n\n', allBV(indBV));
  for bInd= 1:totExp;	      % iteration for different datasets
    % generating data
    [xTrain, yTrain] = sincdata(nType,nTrain,sig02,[],0);
    % initialising the GP
    ogp(1, ...	      	      % input dimension
	1, ...	      	      % output dimension
	'sqexp', ...        % kernel type
	log([1/100 1 1]));    % kernel parameter
    % Assigning a Likelihood.
    % For regression the choices are 'reg_lapl' and 'reg_gauss'.
    ogpinit(@c_reg_gauss,...   % address for log-evid and its ders.
	    1e-4)%,...	      % likelihood (hyper-) parameters
% 	    @em_exp);	      % adjusting likelihood parameters
    net.thresh    = 1e-3;     % the admission threshold for new BVs
    net.maxBV     = maxB;
    net.proj      = 0;

    % HYPERPRIORS used in hyperparameter optimisation.
    ogphypcovpar(1e-2);

    % Initialising GPOPT - optimisation and display options
    gpopt       = defoptions; % NO log-average
    gpopt.ptest = 1;	      % YES test error computation
    gpopt.xtest = xTest;      % the test inputs
    gpopt.ytest = yTest;      % desired outputs
    gpopt.disperr=0;
    gpopt.erraddr={@err_abs}; % function that measures the test error
    gpopt.freq  = nTrain;     % frequency of measuring the errors
    
    gpopt.postopt.isep  = 1;  % USING the TAP/EP algorithm
    gpopt.postopt.itn   = 2;  % number of EP iteration with CHANGING BVs
    gpopt.postopt.fixitn= 1;  % FIXING the BV set.

    % the control of HYPPAR optimisation
    gpopt.covopt.opt    = foptions;% default options
    gpopt.covopt.opt(1) = 1;  % display values
    gpopt.covopt.opt(9) = 1;  % CHECK/NOT gradients
    gpopt.covopt.opt(2) =1e-7;% precision in computing f. value
    gpopt.covopt.opt(3) =1e-7;% display values
    gpopt.covopt.opt(14)=totSCG;% number of iterations
    gpopt.covopt.fnopt ='conjgrad';

    % setting up indices
    lE = 0;
    % training the GP
    for iHyp = 1:totHyp;
      ogptrain(xTrain,yTrain);
      % recomputation of the posterior is NEEDED since OGPTRAIN changed
      % the process parameters.
      ogpreset;
      ogppost(xTrain,yTrain);
      % storing the test errors
      tt = gpopt.testerror; gpopt.testerror = [];
      testErr(indBV,bInd,(lE+1):(lE+length(tt))) = tt;
      lE = lE + length(tt);
      % computing the data evidence
      iTAP(iHyp) = lE;
      eTAP(indBV,bInd,iHyp) = ogpevid([])./nTrain;
      % plotting it
      [meanT, varT] = ogpfwd(xTest);
      stdT          = sqrt(varT);
      meanBV        = ogpfwd(net.BV);
      clf;
      subplot 211; hold on; box on;
      plot(xTest,yTest,'Linewidth',1.5);
      plot(xTest,meanT,'k--','Linewidth',1);
      plot(xTest,meanT + stdT,'r-.','Linewidth',1);
      plot(xTest,meanT - stdT,'r-.','Linewidth',1);
      plot(net.BV,meanBV,'k*','Markersize',12);
      plot(xTrain,yTrain,'r.','Linewidth',0.5);
      title(sprintf('K:%7s, BV#: %d',net.covarfn,size(net.BV,1))); 

      subplot 212; hold on; box on;
      mTAP = permute(mean(eTAP,2),[3 1 2]);
      [ax,h1,h2] = plotyy([1:lE],permute(mean(testErr,2),[3 1 2]),...
	     iTAP,mTAP);
      set(h1,'Linestyle','--','Linewidth',2);
      set(h2,'Linestyle',':','Marker','d','MarkerSize',7);
      legend(ax(1),'Avg. Test Error',2);
      legend(ax(2),'- Avg.Log. Evid.',1);
      title([nType ' likelihood. Lik. par: ' num2str(net.likpar) ]);
      minEv = min(mTAP)-0.01; maxEv = max(mTAP)+.01;      
      ylim(ax(1),[0 .4]); 
      if isfinite(maxEv-minEv);
	ylim(ax(2),[minEv maxEv]);
	set(ax(2),'YTick',linspace(minEv,maxEv,7));
      end;
      xlim(ax(1),[1, lE]);xlim(ax(2),[1, lE]);
      set(ax(2),'XTick',iTAP,'XTickLabel',[1:length(iTAP)]);
      set(ax(1),'XTick',[],'YTick',linspace(0,0.4,7));
      drawnow;
      fprintf('Avg. -ln.ev: %f\n', eTAP(indBV,bInd,iHyp));
      fprintf('Hyper-parameter optimisation #%d\n',iHyp);
      fprintf('K. type:%7s; A: %4.3f, O: %4.3f\n', ...
	      net.covarfn,exp(net.kpar(2:3)));
      fprintf('Input scaling: %4.1f; #BV: %d\n', ...
	      exp(-net.kpar(1)),size(net.BV,1));
      fprintf('Likelihood noise: %10.4f\n\n',net.likpar);
    end;
  end;
end;

