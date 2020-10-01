function [K1, K2] = ogpparadj(K1,K2,mX,sigX2,tolL,tolS)
%OGPPARADJ Adjusts GP such that the TAP/EP it. is not ill-conditioned.
%
%	Description
%
%	[K1, K2] = OGPPARADJ(K1,K2,MX,SIGX2) adjusts the online update
%	parameters such that the resulting GP will not have strongly negative
%	variances.
%
%	The procedure is important in situations where the the likelihood
%	function is not log-concave.  For these cases the  thus the
%	optimisation might not lead to unique solution.  The script changes
%	the update coefficients such that the resulting stochastic process
%	has a positive variance.
%
%	Parameters:
%
%	 K1,K2  - the online update coefficients.
%
%	 MX,SIGX2  - the GP marginals.
%
%	See also
%	OGP, OGPTRAIN
%

%	Copyright (c) Lehel Csato (2001-2004)

% checking the input parameters
if nargin<4;
  error('Not enough parameters');
elseif nargin < 5;
  tolL    = 1e5;
  tolS    = 1e-4;
elseif nargin<6;
  tolS    = 0.00001;
end;

nout     = size(sigX2,2);
[ku,kv]  = eig(-K2); ku = diag(ku);
sK       = find(ku<tolS);
if length(sK)
  fprintf('Expanding K2\n\n');
  ku(sK) = tolS;
  K2     = -kv*diag(ku)*kv';
end;

% adjusting K2 - to make the INFERENCE alg. more stable
[tV,tu] = eig( - (eye(nout)+K2*sigX2)\K2 );
tu      = diag(tu);
iSmall  = find(tu<tolS);
iLarge  = find(tu>tolL);

K2= K2 + eps/100;
c = (K1'/K2)*K1;

if length(iSmall)
%  fprintf('Sml.L.');
  tu(iSmall) = tolS;
end;

if length(iLarge);     % ESSENTIAL FOR S_T_A_B_I_L_I_T_Y
%   fprintf('Lrg.L.');
  tu(iLarge) = tolL;
end;

K2  = - eye(nout)/(tV*diag(1./tu)*tV'+sigX2);
K2  = (K2 + K2')./2;

% the following line is to keep the quadratic part of the log-posterior
% likelihood constant. It is important only when e.g. K2 can be large
% positive number -> replacing it makes a !SMALL! negative one - changing
% COMPLETELY the likelihood contribution.

if K1;
  K1  = K1 * sqrt(abs(c./((K1'/K2)*K1)));
end;
