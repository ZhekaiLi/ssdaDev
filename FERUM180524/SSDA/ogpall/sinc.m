function y=sinc(x)
%SINC	Sin(pi*x)/(pi*x) function. (from MATLAB)
%
%	Description
%
%	SINC(X) returns a matrix whose elements are the SINC of the elements
%	of X i.e.
%	     y = sin(pi*x)/(pi*x)    if x ~= 0
%	       = 1                   if x == 0
%	 where X is an element of the input matrix and Y is the resultant
%	output element. (T. Krauss, 1-14-93 - (c) 1988-1999 The MathWorks)
%
%	See: SQUARE, SIN, COS, CHIRP, DIRIC, GAUSPULS, PULSTRAN, RECTPULS.
%
%	(this implementation is needed when generating toy data)
%
%	See also
%	OGP, OGPTRAIN, C_REG_GAUSS, C_REG_EXP,
%	C_REG_LAPL, DEMOGP_REG
%

%	Copyright (c) Lehel Csato (2001-2004)

y=ones(size(x));
iN=find(x);
y(iN)=sin(pi*x(iN))./(pi*x(iN));
