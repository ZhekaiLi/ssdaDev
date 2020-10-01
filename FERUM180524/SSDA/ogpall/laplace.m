function y=laplace(dims,flipS);
%LAPLACE Sampling from a Laplace distribution.
%
%	Description
%
%	[Y] = LAPLACE(DIMS,FLIPS) - generates a matrix with random variables
%	distributed according to a Laplace distribution with parameter
%	LAMBDA=1.
%
%	Function parameters:
%
%	  DIMS  - the dimensions of the returning matrix Y.     FLIPS  - an
%	indicator to switch between symmetric and   nonsymmetric noises
%	(default is 1).
%
%
%	See also
%	OGP, OGPTRAIN, C_REG_EXP, DEMOGP_REG
%

%	Copyright (c) Lehel Csato (2001-2004)

% argument checking
if nargin==1;
  flipS=1;
end;

y  = rand(prod(dims),1);
y0 = find(y==0);
while length(y0);
  y(y0) = rand(length(y0),1);
  y0    = find(y==0);
end;
y = - log(y);
if flipS;
  y = y .* sign(rand(size(y)) - 0.5);
end;

y = reshape(y,dims);
