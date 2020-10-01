%DEMOGP_MATERN Example to use an external covariance function.
%
%	Description
%
%	DEMOGP_MATERN - script to exemplify the use of an external covariance
%	function.  The function chosen is the Matern covariance function
%	which allows the transition from a rough (once- or twice-
%	differentiable) kernel to the infinitely many times differentiable
%	RBF kernel (implemented in this package as 'SQEXP').
%
%	The demonstration starts from sampling from a realisation of a GP
%	with zero mean and a MATERN covariance of order 1.5. This sample
%	function is corrupted with additive Gaussian noise, then a GP with
%	matern kernel - but different parameters - is trained on the samples.
%
%	In the demo script the experiments are repeated three times and
%	different BV set sizes are also tested (to make it faster, change the
%	parameters accordingly).
%
%	See also
%	OGP, OGPTRAIN, COV_MATERN, SINCDATA
%

%	Copyright (c) Lehel Csato (2001-2004)

% clearing the workspace 
clear all;
global net gpopt ep;

% General information
nTest  = 201;		      % size of the test set
nTrain = 650;		      % size of the training set.
sig02  = 0.01;		      % the noise variance
allBV  = 30:70:170;	      % BV set sizes
totExp = 1;		      % number of experiments for EACH BV set
totHyp = 10;		      % total numbers of Hyperparameter estimation

indBV = 0;
for maxB = allBV;
  indBV = indBV + 1;
  fprintf('\n %4d :\n', allBV(indBV));
  randn('state',2425);        % SETTING the same experiments
  rand( 'state',2245);
  
  for bInd= 1:totExp;	      % iteration for different datasets
    fprintf('\nExperiment #%d\n\n',bInd);
    % initialising the GP
    ogp(1, ...	      	      % input dimension
	1, ...	      	      % output dimension
	'user', ...           % kernel type
	log([1./80 10 1.5]));
    % for user-specified kernels one HAS to define the way to compute the
    net.kfnaddr   = @cov_matern;% KERNEL
    net.gradkaddr = @covgrad_matern;% GRADIENT wrt parameters
    net.thresh    = 1e-4;     % the admission threshold for new BVs
    net.maxBV     = maxB;
    % DEFINING the likelihood
    ogpinit(@c_reg_gauss,...  % address for log-evid and its ders.
	    0.1, ...   	      % likelihood (hyper-) parameters
	    @em_gauss);	      % adjusting likelihood parameters
    % SETTING the distribution of the priors to the parameters:
    % HYPERPARAMETERS
    ogphypcovpar(5e-2);

    % GENERATING THE DATA
    xTest = linspace(-20,20,nTest)';
    xTrain= 40*(rand(nTrain,1)-0.5);
    yTest = ogpsample([xTest; xTrain]);
    yTrain= yTest([(nTest+1):(nTest+nTrain)],:);
    yTest = yTest(1:nTest);
    yTrain = yTrain + sqrt(sig02)*randn(size(yTrain));

    % CHANGING THE INITIAL PARAMETERS OF THE _OGP_
    net.kpar = log([1./5 1 15]);
    % clearing the structure (needed if multiple BV sets are tested)
    ogpreset

    % Initialising GPOPT - optimisation and display options
    gpopt       = defoptions; % NO log-average
    gpopt.ptest = 1;	      % YES test error computation
    gpopt.xtest = xTest;      % the test inputs
    gpopt.ytest = yTest;      % desired outputs
    gpopt.disperr=0;
    gpopt.erraddr={@err_mse};   % function that measures the test error
    gpopt.freq  = floor(nTrain);% frequency of error meas.
    
    gpopt.postopt.isep  = 1;  % USING the TAP/EP algorithm
    gpopt.postopt.itn   = 3;  % number of EP iteration with CHANGING BVs
    gpopt.postopt.fixitn= 1;  % FIXING the BV set.

    gpopt.covopt.opt   = foptions;
    % setting the FOPTIONS indicator to check the gradients
    gpopt.covopt.opt(1)  =  -1;% DISPLAY/NOT training values
    gpopt.covopt.opt(9)  =  0;% CHECK/NOT   gradients
    gpopt.covopt.opt(2)   = 1e-3;% precision in computing optimum values
    gpopt.covopt.opt(3)   = 1e-3;
    gpopt.covopt.opt(14)  = 20;
    gpopt.covopt.fnopt = 'scg';

    % saving the CURVES themselves for later reference
    saveYTest(indBV,bInd,:) = yTest;
    
    % setting up indices
    lE = 0;
    % training the GP
    for iHyp = 1:totHyp;
      ogptrain(xTrain,yTrain);
      saveBV(indBV,bInd,iHyp) = size(net.BV,1);
      ogpreset;
      % saving the OGP configuration
      saveGP(indBV,bInd,iHyp) = net;

      % setting the parameters to a "reasonable value"
      hp = exp(net.kpar);
      if hp(2)<0.3; hp(2) = 0.3; end;
      if hp(2)>200; hp(2) = 200; end;
      if hp(3)<0.7; hp(3) = 0.7; end;
      if hp(3)>15;  hp(3) = 15;  end;
      net.kpar = log(hp);
      if net.likpar>20; net.likpar = 20; end;

      % computing again the posterior approximation
      ogppost(xTrain,yTrain);
      % storing the test errors
      tt = gpopt.testerror; gpopt.testerror = [];
      testErr(indBV,bInd,(lE+1):(lE+length(tt))) = tt;
      lE = lE + length(tt);
      % plotting it
      [meanT, varT] = ogpfwd(xTest);
      stdT          = sqrt(varT);
      meanBV        = ogpfwd(net.BV);
      clf;
      axes('position',[0.1 0.35 0.82 0.57]); hold on; box on;
      plot(xTest,yTest,'Linewidth',1.5); plot(xTest,meanT,'k--','Linewidth',1);
      plot(xTest,meanT + stdT,'r-.','Linewidth',1);
      plot(xTest,meanT - stdT,'r-.','Linewidth',1);
      plot(net.BV,meanBV,'k*','Markersize',12);
      plot(xTrain,yTrain,'r.','Linewidth',0.5);
      
      title(sprintf('X sc.: %3.2f; K: (%3.2f, %3.2f); L. var: %4.3f', ...
		    exp(-net.kpar(1)),exp(net.kpar(2:3)), ...
		    net.likpar));
      set(gca,'XTick',[]);

      axes('position',[0.1 0.1 0.82 0.19]); hold on; box on;
      semilogy(permute(sum(testErr,2),[3 1 2]));
      ylabel('Test Error');
      legend(num2str(allBV'));
      ylim([0 0.03]);
      set(gca,'Ygrid','on','YTick',[0 0.002 0.005 0.01 0.2]);
      drawnow;

      % computing the data evidence
      eTAP = ogpevid([]);
      fprintf('evidence: %f;  BV#: %3d\n',eTAP,size(net.BV,1));
      fprintf('X scale: %3.2f; ',exp(-net.kpar(1)));
      fprintf('K(A: %3.2f, O: %3.2f)  ',exp(net.kpar(2:3)));
      fprintf('L.var: %4.3f\n\n',net.likpar);
    end;
  end;
end;

% ONE CAN display the statistics/evolution of different hyperparameters
% from the saveGP structure
save results_matern
