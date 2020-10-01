function Evid = ogpevid(hp)
%OGPEVID Evaluates the evidence for Sparse OGP.
%
%	Description
%
%	EVID = OGPEVID with the use of the GLOBAL Gaussian Process data
%	structure NET and -- if present -- the parameters of the TAP/EP
%	structure EP, the function returns the negative log-evidence of the
%	data.
%
%	Prior to evaluation of the evidence the parameters of the model are
%	recomputed.  The update does not occur if HP is empty or there is no
%	input specified.  The format of the parameter vector hp is specified
%	by OGPPAK.
%
%	If the structure EP is not empty, then the computation of the
%	evidence includes terms from the TAP/EP approximation scheme.
%
%	If the TAP/EP representation of the Gaussian process is not given,
%	then the returned value is only approximate.
%
%	See also
%	OGP, OGPCOVARF, OGPPAK, OGPEVIDGRAD
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% update NET
if nargin>0 & ~isempty(hp);  ogpunpak(hp); end;

nBV  = size(net.BV,1);
if ~nBV;		      % LARGE NUMBER
  Evid = 1e10; return;
end;

facM = exp(net.kpar(2)-50+log(net.addEps));
KB   = net.KB+eye(nBV)*facM;

if ~isfield(ep,'lamP');% no EP coefficients are available
  % fprintf('NO TAP/EP coefficients! Evidence is approximate\n');
  Z = eye(nBV) + KB *net.C;
  Evid = sum(log(eig(Z,'nobalance'))) ...
	 + (net.w'/Z)*(KB*net.w);
else
  nZ   = full(diag(ep.lamP));
  nSel = nZ(find(nZ>eps/100));
  Evid = sum(log(nSel)) - sum(ep.aP.^2.*nZ) + ...
	 2*sum(ep.logZ) - length(nSel)*log(2*pi);
  Ksm  = (ep.projP'*ep.lamP*ep.projP)*KB + eye(nBV);
  kall = ep.projP'*ep.lamP *ep.aP;
  Evid = Evid + kall'*(KB/Ksm)*kall;
  Evid = Evid - sum(log(eig(Ksm,'nobalance')));
end;

Evid = - Evid/2;

% implementing priors for the hyperparameters
if find(net.hyplambda);
  pars = (ogppak - net.hypmean).^2;
  Evid = Evid + net.hyplambda*pars'/2;
end;

% $$ NET.ADDEPS - defines an extra "ridge" to make the kernel matrix
% pos.def. => one could ask for putting a penalty for non-zero addEps
% values -> THIS IS NOT NEEDED - OGP WILL REDUCE THE NUMBER OF BVs.
% $$ NET.W and NET.C are all using an extra inversion -- using them to
% compute the evidence IS NOT ADVISED.
