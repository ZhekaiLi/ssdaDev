function ogpstep_sp(indX,kX,K1,K2, gamma, hatE);
%OGPSTEP_SP Performs a sparse online update step of the Gaussian Process.
%
%	Description
%
%	OGPSTEP_SP(INDX,KX,K1,K2,GAMMA,HATE) - updates the GLOBAL Gaussian
%	process structure NET based on information about the current input.
%	This function is called from OGPPOST.
%
%	In the sparse update step the set of basis vectors is not changed!
%	The current input modifies the only the coefficients of the GP.
%
%	See the online (HTML) documentation for details
%
%	See also
%	OGP, OGPTRAIN, OGPSTEP_FULL
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% computing scaling factor \eta
if net.proj;		      % PROJ. update: keeps can.par
  eta = eye(net.nout) + K2*gamma;
else			      % MARG. update: keeps marg.
  eta = eye(net.nout);
end;

% updating GP parameters
tVector = net.C*kX + hatE;
net.w   = net.w + tVector*(eta\K1);
eta     = eta\K2;
net.C   = net.C + tVector*eta*tVector';
net.C   = (net.C + net.C')./2;

% UPDATE of the ep-parameters is made in !ogpstep_ep!
