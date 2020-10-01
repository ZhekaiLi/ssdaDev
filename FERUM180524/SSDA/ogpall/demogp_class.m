%DEMOGP_CLASS Two-dimensional binary classification example.
%
%	Description
%
%	DEMOGP_CLASS - demo script for using the OGP package for
%	classification.
%
%	The script generates data from a mixture of three Gaussians with two
%	of them belonging to 'CLASS1' and the third to 'CLASS2'.
%
%	After the generation of the training data, iterates the computation
%	of the posterior (with the sparse online learning algorithm) and the
%	adjustment of the kernel hyperparameters in the OGP (net).
%
%	Similar to regression (DEMOGP_REG), the script computes and displays
%	average test errors for different BV set sizes.
%
%	The kernel family and the number of iteratins have to be set manually
%	in the script - see DEMOGP_CLASS_GUI for a more interactive
%	classification demo.
%
%	See also
%	OGP, OGPTRAIN, C_CLASS_BIN, DEMOGP_CLASS_GUI, DEMOGP_REG
%

%	Copyright (c) Lehel Csato (2001-2004)

% DEMONSTRATION OF THE ONLINE GAUSSIAN PROCESS PACKAGE FOR CLASSIFICATION
% 
% This matlab script illustrates the usage of the OGP package for binary
% classification. 

clear all; close all; clc
global net gpopt ep;

% other data
allBV    = 50:20:50;	      % upper limit for the basis set.
totHyp   = 8;		      % number of iterations for Hyperp.opt.
totSCG   = 6;

% datasets
nTest       = 500;
nTrain      = 150;
% generating test data on a grid.
gSize = ceil(sqrt(nTest));
[xTest,yTest] = cl_data(gSize,[],1);
nTest = gSize^2;

indBV = 0;
for maxB = allBV;
  indBV = indBV + 1;
  fprintf('Classification using  %4d BVs:\n', allBV(indBV));
  
  % generating data
  [xTrain, yTrain] = cl_data(nTrain);
  % building a GP
  ogp(2, ...	      	      % input dimension
      1, ...	      	      % output dimension
      'poly', ...          % kernel type
      log([1/10 1/10 6 6]));  % kernel parameter
  ogpinit(@c_class_bin,5);
  net.thresh    = 2e-3;     % the admission threshold for new BVs
  net.maxBV     = maxB;
  % HYPERPRIORS used in hyperparameter optimisation.
  ogphypcovpar(1e-1);
  net.hyplambda(end) = 0.001;
  net.hypmean(end) = net.kpar(end);

  % Initialising GPOPT - optimisation and display options
  gpopt               = defoptions;
  gpopt.postopt.isep  = 1;  % USING the TAP/EP algorithm
  gpopt.postopt.itn   = 2;  % number of EP iteration with CHANGING BVs
  gpopt.postopt.fixitn= 2;  % FIXING the BV set.
  
  gpopt.covopt.opt    = foptions;% default options
  gpopt.covopt.opt(1) = -1;   % display values
  gpopt.covopt.opt(9) =  1;   % CHECK/NOT gradients
  gpopt.covopt.opt(2) =1e-7;  % precision in computing f. value
  gpopt.covopt.opt(3) =1e-7;  % display values
  gpopt.covopt.opt(14)=totSCG;% number of iterations
  gpopt.covopt.fnopt  ='conjgrad';
  clf;
  % training the GP
  for iHyp = 1:totHyp;
    %keyboard
    ogptrain(xTrain,yTrain);
    ogpreset;
    ogppost(xTrain,yTrain);
    savePars(:,iHyp,indBV) = net.kpar(1:end-1);
    saveEV(iHyp,indBV)     = ogpevid([]);

    hypSet = exp(ogppak);
    hypset(2:3) = sqrt(1./hypSet(2:3));
    fprintf(['\n\nIt #%d (%d BV) parameters:\nBias: %2.1f\nInput ' ...
	     'scale: (%4.3f, %4.3f)\nKernel amp.: %4.3f,  '...
	     'ord.: %2.1f\n\n'],iHyp,size(net.BV,1),hypSet);
    % measuring the errors
    [meanT, varT] = ogpfwd(xTest);
    stdT          = sqrt(varT);
    predP = (1+erf(meanT./(sqrt(2)*stdT)))/2;
    cla; hold on; box on;
    % plotting the data from Class 1
    iInd = find(yTrain==1);
    hTrainP = plot(xTrain(iInd,1),xTrain(iInd,2),'mv');
    iInd = find(yTrain==-1);
    hTrainN = plot(xTrain(iInd,1),xTrain(iInd,2),'c+');
    % plotting the BVs
    hBV = plot(net.BV(:,1),net.BV(:,2),'kd','Linewidth',3,'MarkerSize',12');
    tX1 = reshape(xTest(:,1),gSize,gSize);
    tX2 = reshape(xTest(:,2),gSize,gSize);
    tY  = reshape(yTest     ,gSize,gSize);
    tG  = reshape(predP     ,gSize,gSize);
    tV  = reshape(stdT      ,gSize,gSize);
    % Class boundaries -  Theoretical - the one obtained by the GP
    if strcmp(version('-release'),'14');
      [ttt, hTh] = contour('v6',tX1(1,:),tX2(:,1),tY,[.5,.5],'k--');
      [ttt, hGP] = contour('v6',tX1(1,:),tX2(:,1),tG,[.5 .5],'r');
    else
      [ttt, hTh] = contour(tX1(1,:),tX2(:,1),tY,[.5,.5],'k--');
      [ttt, hGP] = contour(tX1(1,:),tX2(:,1),tG,[.5 .5],'r');      
    end;
    axis([-2.5 2.5 -2.5 2.5]);
    title(['BV set size ' num2str(size(net.BV,1))]);
    legend([hTrainP, hTrainN, hBV, hTh, hGP'], ...
	   'Data Class 1','Data Class 2','Basis Vectors', ...
	   'Theory','GP pred');
    drawnow;    
  end;
end;
