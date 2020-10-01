%DEMOGP_FIXED Gaussian Process regression using a fixed set of basis vectors.
%
%	Description
%	DEMOGP_FIXED - script that implements Online Gaussian Process
%	regression with a fixed set of BV set.   The script generates data
%	from the noisy SINC function where the  inputs are 2-dimensional and
%	the radius of the input is the argument to SINC.
%
%	In this script the BV set is taken from the test set, thus
%	independent  from the training set and no addition or removal to/from
%	the BV set  is made during learning - using the sparse online
%	algorithm.
%
%	The script performs hyper-parameter optimisation and the performance
%	of the fixed-size OGP is compared with the "normal" algorithm using
%	different BV set sizes.
%
%	Similar to one-dimensional regression, the script computes test error
%	and displays the averages over different runs. Additionally, it
%	displays the  position of the BV elements and the posterior mean. If
%	visual inspection of the results is wanted, then one must set the
%	variable ISBREAK to a nonzero value.
%
%	See also
%	OGP, OGPTRAIN, DEMOGP_REG
%

%	Copyright (c) Lehel Csato (2001-2004)

% INITIALISING THE WORKSPACE
clear all; clc;
global net gpopt ep;

% FOR an interactive demonstration, set isBreak to a nonzero value
isBreak = 1;

fprintf(['Sparse Gaussian Process regression with fixed BV set.\n' ...
	 'The BV set is the TEST set, as in ' ...
	 'TRANSDUCTIVE LEARNING.\n\nThe initial ' ...
	 'set of BVs is changed as a result of\n' ...
	 'hyper-parameter optimisation: the inefficient\nBVs ' ...
	 'are removed\n\n' ]);
fprintf(['In the first part the OGP structure is initialised to the ' ...
	 'TEST set\nThat is, all OGP parameters are zero except for ' ...
	 'their LOCATION.\nThe inference and model selection ' ...
	 'steps alternate to obtain\nthe optimal Gaussian process.\n']);
fprintf(['At the end of each iteration the current hyperparameters are\n' ...
	 'printed out and the figure shows the true/appr. functions.\n\n']);
fprintf(['A second step uses the conventional BV set selection from ' ...
	 'the\ntraining set and the hyperparameters from the first ' ...
	 'part and the\ntest errors for different BV set sizes are ' ... 
	 'computed.\n\n']);

if isBreak; 
  fprintf('Press any key to continue.\n'); pause; 
else 
  fprintf('To pause the program, set ''isBreak'' in the script.\n'); 
end;

%____________________________________________________________ 
% SETTING up the different parameters for test/training data
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sig02 = 0;		          % noise variance
nType = 'gauss';	      % noise TYPE
nTest   = 100;		      % TEST set size == BV SET size
nTrain  = 1000;		      % TRAINING set size
nSqrt   = sqrt(nTest);

allBV    = 50:10:50;	      % BV set sizes for comparison of the
                              % transductive learning with normal case
totExp   = 1;		      % number of experiments
totHyp   = 5;

%____________________________________________________________ 
% INITIALISING THE GP
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ogp(2,1, ...	      	      % input & output dimension
    'sqexp', ...       	      % kernel type
    log([.1 .10 ...   	      % length-scales
	 1 1]));      	          % amp. & order
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ASSIGNING A LIKELIHOOD FUNCTION
ogpinit(@c_reg_gauss,...      % address of the likelihood function
	1e-4 );
net.thresh    = 1e-4;        % the admission threshold for new BVs
net.maxBV     = 1000;
ogphypcovpar(5e-2);
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% SAVING THE ORIGINAL PARAMETERS
netS = net;

% SETTING THE TEST POINTS
[xTest,  yTest ] = sinc2data(nType,nTest,0,[],1,20,20);

% Initialising GPOPT - optimisation and display options
gpopt       = defoptions;   % NO log-average
gpopt.ptest = 1;	      % YES test error computation
gpopt.xtest = xTest;          % the test inputs
gpopt.ytest = yTest;          % desired outputs
gpopt.disperr=0;
gpopt.erraddr={@err_mse};     % function that measures the test error
gpopt.freq  = 10;             % frequency of measuring the errors
gpopt.postopt.isep   = 1;     %  USING the TAP/EP algorithm
gpopt.postopt.itn    = 3;     % number of EP iteration with CHANGING BVs
gpopt.postopt.fixitn = 1;     % FIXING the BV set.
gpopt.covopt.opt     = 2;
gpopt.covopt.fnopt   = 'scg';
gpopt.postopt.isrm   = 1;     % flag to remove BV or not
% SAVING OPTIONS
gpoptS = gpopt;

%____________________________________________________________ 
% PARAMETERS for visualisation
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gX       = [-20:0.5:20];
gN       = length(gX);
[pX pY]  = meshgrid(gX);
pX       = [pX(:) pY(:)];

for bInd= 1:totExp;	      % iteration for different datasets
  % generating data
  % training points are random (encoded with 0 as the last parameter)
  [xTrain, yTrain] = sinc2data(nType,nTrain,sig02,[],0,20,20);

  net   = netS;
  gpopt = gpoptS;
  
  % FIXED BV set - INITIALISATION
  ogpemptybv(xTest(randperm(size(xTest,1)),:));
  netF.isBVfixed = 1;

  % perform training for fixed BV set
  fprintf('\nExperiment #%d\n',bInd);
  for iHyp=1:totHyp;
      if iHyp > 5; gpopt.postopt.isrm   = 0; end
    ogptrain(xTrain,yTrain);
    ogpreset;  % remove the BV
%     % obtaining the posterior
    ogppost(xTrain,yTrain);
    disp('mean error'); disp(mean(gpopt.testerror));
    fixErr(bInd,:) = gpopt.testerror; gpopt.testerror = [];
    
    % plotting the mean of the GP using the fixed size GP
    figure(1); clf; hold on;
    pM = ogpfwd(pX);
    subplot('position',[.05 .05 0.55 0.85]); hold on;
    surf(gX,gX,reshape(pM,gN,gN),'EdgeAlpha',0.3);
    % fancy plot of the original surface
    x1 = reshape(xTest(:,1),[nSqrt nSqrt]);
    x2 = reshape(xTest(:,2),[nSqrt nSqrt]);
    y  = reshape(yTest,    [nSqrt nSqrt]);
    h  = surf(x1,x2,y,0.5*ones(size(y))); view(3);
    set(h,'FaceAlpha',0.3,'EdgeAlpha',0.5,'Linewidth',0.5);
    title('Regression with fixed BVs');
    zlim([-5 20]); grid on;
    subplot('position',[.62 .02 .33 .88]);
    h=scatter(net.BV(:,1),net.BV(:,2),'rp');
    axis equal;
%     set(h,'Markersize',14,'Linewidth',2);
    hold on; box on; axis([-20 20 -20 20]);
    title(['BV positions (#BV=' num2str(size(net.BV,1)) ')']);
    h=scatter(xTrain(:,1),xTrain(:,2),'kx');
%     set(h,'Markersize',10,'Linewidth',1);
    drawnow;
    % computing the data evidence
    eTAP = ogpevid([]);
    fprintf('evidence: %f\n',eTAP);
    fprintf('Current hyper-parameters:\n');
    fprintf('Input scaling: %4.1f %4.1f ; BV set size: %d\n', ...
	    exp(-net.kpar(1:2)),size(net.BV,1));
    fprintf('Kernel amplitude: %7.4f\n',exp(net.kpar(1)));
    fprintf('Likelihood noise: %10.4f\n\n',net.likpar);
    if isBreak; fprintf('Press any key to continue.\n'); pause; end;
  end;
  % SAVING the process/parameters
  netF = net;
  
  fprintf(['In the second part of the demo, using the hyper-parameter\n' ...
	   'set found, we train the OGP using a fixed BV set\n\n']);
  if isBreak; fprintf('Press any key to continue.\n'); pause; end;
  
  % perform training for varying BV set sizes
  indBV = 0;
  for maxB = allBV;
    indBV = indBV + 1;
    fprintf('\n Training with %4d BV...\n', allBV(indBV));
    % training the GP
    net        = netS;
    net.kpar   = netF.kpar;   % copy over hyperparameters
    net.likpar = netF.likpar;    
    net.maxBV  = maxB;
    gpopt      = gpoptS;

    for iHyp = 1:totHyp;
      ogpreset;
      ogptrain(xTrain,yTrain);
      ogpreset;
      ogppost(xTrain,yTrain);
      chErr(indBV,bInd,:) = gpopt.testerror; gpopt.testerror = [];
    end;
    % computing the data likelihood
    eTAP = ogpevid([]);
    fprintf('evid: %f\n',eTAP);
    % plotting cumulative statistics
    figure(2); clf; hold on; box on;
    plot(mean(fixErr,1),'k--');
    plot(permute(mean(chErr,2),[3 1 2]));
    ylim([0 5]);
    legend('Test error using a fixed set', ... 
	   [num2str(allBV(1)) ' Basis Vectors'], ...
	   [num2str(allBV(2)) ' Basis Vectors'], ...
	   [num2str(allBV(3)) ' Basis Vectors'])
    drawnow;
  end;

  if isBreak; fprintf('Press any key to continue.\n'); pause; end;

end;
