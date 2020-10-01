%DEMOGP_SIMP Gaussian process regression using user inputs.
%
%	Description
%
%	DEMOGP_SIMP - is an interactive demo of the Gaussian process
%	regression.
%
%	The code is simple, designed to provide a visual demonstration tool
%	where the Gaussian Process learns the model as new data is presented
%	to it.
%
%	Two windows show up: the left one is for adding training data and in
%	the right window samples from the posterior are drawn.
%
%	The left figure shows the mean of the posterior process, two standard
%	deviation of the posterior and two standard deviation of the
%	predictive distribution (the latter one is displayed assuming
%	Gaussian likelihood model).
%
%	Data addition is ended if return is pressed instead a mouse click
%	when adding an input.
%
%	The parameters of the Gaussian process - the type of kernel and the
%	likelihood function - have to be set in the source code.
%
%	See also
%	OGP, OGPTRAIN, DEMOGP_REG_GUI
%

%	Copyright (c) Lehel Csato (2001-2004)


% clearing the workspace
clear all;
global net gpopt ep;

% initialising the GP
ogp(1, ...	      % input dimension
    1, ...	      % output dimension
    'ratquad' , ...	      % kernel type
    log([1/15 1 6]));   % kernel parameter
% HYPERPRIORS used in hyperparameter optimisation.
ogphypcovpar(5e-3);

% assigning a likelihood to the data
ogpinit(@c_reg_gauss,...      % address for log-evid and its ders.
	0.025, ...            % likelihood (hyper-) parameters
	@em_gauss);
net.thresh    = 2e-3;         % the admission threshold for new BVs
net.maxBV     = 60;
% experiments using the BIAS
net.bias         =  -10;
net.hyplambda(1) = 1e-8;
net.hypmean(1)   =  -10;

% Initialising GPOPT - optimisation and display options
gpopt       = defoptions; % NO log-average
gpopt.postopt.isep  = 1;      % USING the TAP/EP algorithm
gpopt.postopt.itn   = 2;      % number of EP iteration with CHANGING
gpopt.postopt.fixitn= 1;      % FIXING the BV set.

% parameters to the (NETLAB) optimisation algorithm
gpopt.covopt.opt    = foptions;
gpopt.covopt.opt(1) =  -1;    % DISPLAY/NOT training values
gpopt.covopt.opt(9) =   0;    % CHECK/NOT   gradients
gpopt.covopt.opt(2) =1e-3;    % precision in computing optimum values
gpopt.covopt.opt(3) =1e-3;
gpopt.covopt.opt(14)=   floor(1.5*length(ogppak));    % # of steps
gpopt.covopt.fnopt  = 'conjgrad';

% setting up the plots
xTrain = []; yTrain =[];
xTest  = [-5:0.1:5]';
h_1 = figure(1);
set(h_1,'units','normalized');
set(h_1,'position',[0 0.25 0.5 0.7]);
clf; axes('position',[0.03 0.03 0.94 0.90]); box on;
h_2 = figure(2); 
set(h_2,'units','normalized');
set(h_2,'position',[0.5 0.25 0.5 0.7])
clf; axes('position',[0.03 0.03 0.94 0.94]); box on;
axis([-5 5 -5 5]);

% performing iterations
while 1;
  % generating data
  figure(h_1);  axis([-5 5 -5 5]);
  tp = ginput(1);
  if isempty(tp); return; end;
  xTrain = [xTrain; tp(1)];
  yTrain = [yTrain; tp(2)];
  % resetting the posterior is NEEDED
  ogpreset;
  % more than a single iteration to fit the data
  for iFit = 1:5;
    ogptrain(xTrain,yTrain);
    ogpreset;
  end;
  % computing the posterior
  ogppost(xTrain,yTrain);
  fprintf(['%2d BV, bias: %-4.1f, K: %s (a: %4.2f,o: %4.2f), ' ...
	   'I.sc: %4.1f; L.d: %5.3f, -LnEv: %5.2f\n'],...
	  length(net.BV),net.bias,...
	  net.covarfn,exp(net.kpar(2:3)), exp(-net.kpar(1)), ...
	  sqrt(net.likpar),ogpevid([])/length(xTrain));
  % plotting the posterior
  [meanT, varT] = ogpfwd(xTest);
  stdT          = sqrt(varT);
  stdPr         = sqrt(varT + net.likpar);
  meanBV        = ogpfwd(net.BV);
  figure(h_1);
  xlim([-5 5]); ylim([-5 5]);
  cla;  box on; hold on;
  plot(xTest,meanT,'b--','Linewidth',4);
  plot(xTest,meanT + 2 * stdT,'r-.','Linewidth',2.5);
  plot(xTest,meanT - 2 * stdT,'r-.','Linewidth',2.5);
  plot(xTest,meanT + 2 * stdPr,'g-','Linewidth',2);
  plot(xTest,meanT - 2 * stdPr,'g-','Linewidth',2);
  plot(xTrain,yTrain,'k+','Linewidth',2.5,'Markersize',18);
  plot(net.BV,meanBV,'y*','Markersize',12,'Linewidth',3);
  title(sprintf('Tr#: %3d, BV#: %2d, L.dev: %4.3f',length(xTrain), ...
		length(net.BV),sqrt(net.likpar)),'FontSize',16);
  figure(h_2); cla;
  xlim([-5 5]); ylim([-5 5]); box on; hold on;
  % plotting the random samples from OGP
  for ii=1:15;
    yT(:,ii) = ogpsample(xTest);
  end;
  plot(xTest,yT,'-','Linewidth',2);
  drawnow;
  % SAVING THE DATA TO COMPARE WITH BATCH
    save dataP xTrain yTrain net gpopt
end;

