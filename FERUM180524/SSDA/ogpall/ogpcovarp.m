function [covp, covf] = ogpcovarp(x1, x2)
%OGPCOVARP Calculate the prior covariance for the Sparse GP.
%
%	Description
%
%	[COVP, COVF] = OGPCOVARP(NET, X1, X2) considers the global OGP data
%	structure NET and the two matrices X1, X2 of input vectors and
%	computes the matrix of the prior covariance.  This is the function
%	component of the covariance plus the exponential of the bias term. If
%	called with a single input X1, the function returns only the diagonal
%	of OGPCOVARP(X1,X1).
%
%	[COVP, COVF] = OGPCOVARP(X1, X2) also returns the function component
%	of the covariance.
%
%	Relation between the function component and the value returned by
%	OGPCOVARP is
%	covp = exp(net.bias) + covf
%	 where COVF is the covariance matrix returned by OGPCOVARF and bias
%	is added to each kernel element.
%
%	For the available covariance functions see OGPCOVARF.
%
%	The variable NET is global: the kernel parameters and the bias value
%	are taken from this structure.
%
%	See also
%	OGP, OGPINIT, OGPCOVARF
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% consistency checking
errstring = consist(net, 'ogp', x1);
if ~isempty(errstring);
  error(errstring);
end

if nargin==1;		      % only computes the diagonal covariance
  covf  = ogpcovarf(x1);
else
  if size(x1, 2) ~= size(x2, 2)
    error('Number of variables in x1 and x2 must be the same');
  end;
  covf = ogpcovarf(x1, x2);
end

covp = covf + exp(net.bias);
