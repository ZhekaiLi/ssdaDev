function [kl] = ogpkl(net2,sameBV);
%OGPKL	Computes the KL-distance of the GP marginals.
%
%	Description
%
%	[KL]=OGPADJGP(NET2) computes the KL-distances between two Gaussian
%	processes which use the same covariance function.  One of the
%	Gaussian processes is given in the GLOBAL structure NET, the second
%	is ginven as a parameter.
%
%	The same kernel function for the two processes is essential:  this is
%	needed to compute the KL-distance between GPs (exact computation is
%	only possible in this case).
%
%	The computation is much faster if the BV sets for the GLOBAL NET and
%	NET2 are the same (i.e. when using a GP inference with fixed BV set
%	or transductive learning). This is indicated with a third boolean
%	argument SAMEBV set to a nonzero value.
%
%	This function can be used as a Cauchy-termination criterion to the
%	online iterations, if the successive OGP structures are close in the
%	KL-sense, this means that the aflgorithm has reached a stable
%	solution.
%
%	See also
%	OGP, OGPTRAIN, OGPCOVARF, OGPCOVARP
%

%	Copyright (c) Lehel Csato (2001-2004)

% GLOBAL STRUCTURE NET is one of the GP-s
global net gpopt ep;

% error checking
if sum(abs(ogppak(net)-ogppak(net2)))>1e-14;
  error(['Cannot compute KL-distance of processes with different' ...
	 ' kernels']);
end;

% input argument checking
if nargin<2;
  sameBV = 0;
end;

% BV set sizes
n1  = size(net.w,1);
n2  = size(net2.w,1);

if sameBV;		      % SIMPLE case: same BVs
  if n1~=n2;
    error('BV''s are not the same for KL-distance');
  end;
  S_s = eye(n2) + net2.C*net2.KB;
  m_d = net2.w - net.w;
  kl  = -log(det(eye(n2) + net.C*net2.KB)) + log(det(S_s));
  S_s = net2.KB/S_s;
  kl  = kl + trace((net.C-net2.C)*S_s) ...
	+ (m_d'*S_s*m_d);
  
else			      % DIFFERENT BVs

  % auxiliary quantities
  S_s = eye(n2) + net2.C*net2.KB;

  % temporary variable - first the LINK matrix
  % since globals are used, NET has to be SAVED ....
  netSaved = net; net      = net2;
  
  Kpr = ogpcovarp(net2.BV,netSaved.BV);
  % ... !THEN! put back
  net = netSaved;

  % second - the distance in the input space
  m_d = net2.KB*net2.w - Kpr*net.w;

  kl  = trace(net.C*net.KB) ...
	- trace(((net2.KB + Kpr*net.C*Kpr')/S_s)*net2.C) ...
	- log(det(eye(n1) + net.C*net.KB)) ...
	+ log(det(S_s)) ...
	+ m_d'*(((net2.KB*S_s)\m_d)) ...
	+ net.w'*(net.KB-Kpr'*(net2.KB\Kpr))*net.w;
end;

kl  = kl/2;
