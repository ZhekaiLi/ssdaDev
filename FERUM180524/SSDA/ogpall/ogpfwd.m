function [y, sigsq] = ogpfwd(x, ktest)
%OGPFWD	Forward propagation through a Sparse OGP.
%
%	Description
%
%	Y_MEAN = OGPFWD(X) using the GLOBAL sparse OGP data structure NET
%	together with the matrix X of input vectors, the function forward
%	propagates the inputs through the model to generate a matrix Y_MEAN
%	of posterior means.  Each row of X corresponds to one input vector
%	and each row of Y_MEAN corresponds to one output vector, the mean of
%	the approximated posterior GP.
%
%	[Y_MEAN, Y_VAR] = OGPFWD(X) also generates a column vector Y_VAR of
%	conditional variances, each row corresponds to the variance of an
%	input pattern.
%
%	If the network output has more than a single dimension, then the
%	outputs (both Y_MEAN and Y_VAR) will be [NOUT*NX,1] and
%	[NOUT*NX,NOUT] respectively.
%
%	[Y_MEAN, Y_VAR] = OGPFWD(NET, X, KTEST) uses a pre-computed
%	covariance matrix KTEST between the inputs X and the BV set in the
%	forward propagation. It assumes that the inputs and the BV set
%	elements do not change. This increases efficiency when several calls
%	to OGPFWD are made.
%
%	See also
%	OGP, OGPINIT, DEMOGP_REG
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% consistency checking
errstring = consist(net, 'ogp', x);
if ~isempty(errstring);
  error(errstring);
end

if nargin<2
  ktest = ogpcovarp(x, net.BV);
end;

% Predict mean
y = ktest*net.w;
if ~isnumeric(net.prmean);
  y = y + feval(net.prmean,x,net.prmeanp);
end;
  
if nargout >= 2   	      % Predict squared variance
  sigsq = ogpcovarp(x);
  if net.nout==1;	      % scalar outputs are simpler
    tt    = sum((ktest*net.C).*ktest,2);
  else
    nX    = size(x,1);
    for iX=1:nX;	      % comp. if we want FULL cov. structure
      tI  = [(net.nout*(iX-1)+1):(net.nout*iX)];
      tkx = ktest(tI,:);
      tt(tI,:) = tkx*net.C*tkx';
    end;
  end;
  sigsq = sigsq + tt;
end
