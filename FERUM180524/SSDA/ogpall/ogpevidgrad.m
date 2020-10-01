function g = ogpevidgrad(hp)
%OGPEVIDGRAD Computes the gradient for the Sparse OGP.
%
%	Description
%
%	G = OGPEVIDGRAD(HP) takes the OGP parameter vector HP and the GLOBAL
%	Sparse OGP data structure NET and returns the gradient of the
%	marginal likelihood (i.e. the gradient of the evidence).
%
%	The gradient is computed only with respect to the kernel parameters,
%	the gradients of the evidence w.r.t. the likelihood parameters are
%	treated differently.  G is a column vector with the dimension of
%	model hyperparameters where the order is provided by OGPPAK.
%
%	This function uses OGPCOVGRAD to compute the gradient values.  If the
%	structure EP is empty, then the gradient (and the evidendce) is only
%	approximate, this might lead to inaccuracies when optimising the
%	kernel parameters.
%
%	The computation of the gradient is can be extremely time and memory-
%	consuming. To alleviate the memory requirements, if there are too
%	many elements to compute, OGPEVIDGRAD makes a sequence of gradient
%	computations, i.e. calls OGPCOVGRAD with an optional argument -- see
%	the documentation of OGPCOVGRAD for details.
%
%	See also
%	OGPCOVARF, OGPPAK, OGPCOVGRAD, OGPTRAIN
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% updates NET
ogpunpak(hp);

isFull = ~isempty(ep);

nbv   = size(net.BV,1);
nbv2  = nbv^2;
ncovp = length(ogppak);
nAll  = nbv2;

% THE gardient wrt. the stored GP only.
gcovBV = net.KB(:) - exp(net.bias);
eG     = reshape(net.w*net.w' + net.C,[nbv2 1])/2;

warning off;
% THE gradient CHANGES if EP is there
if isFull;
  nAll  = size(net.BV,1)*size(ep.X,1);
  if ~isfield(ep,'Kplus');
    ep.Kplus = ogpcovarp(ep.X,net.BV);
  end;
  gcovX = ep.Kplus(:)-exp(net.bias);
  facM   = exp(net.kpar(2)-50+log(net.addEps));
  KB     = net.KB + eye(nbv)*facM;
  KLam   = full(ep.Kplus' * ep.lamP);
  CC     = (KB + KLam * ep.Kplus)\KLam;
  epG    = CC - CC*ep.aP*ep.aP'*ep.lamP + net.w*net.w'*KLam;
  epG = reshape(epG',[nAll,1]);
end;
warning on;

% ONE SHOULD ASK FOR A FULL GRADIENT FOR SMALL TRANING/BV SETS
if nAll*ncovp<1e8;	      % LIMIT the expensive REPMAT operations
  % DERs. wrt. the bias (from ogpcovarp net.bias):
  allG(:,1) = exp(net.bias)*ones(nbv2,1);
  % DERs. wrt. the kernel parameters
  gcovBV          = ogpcovgrad([],gcovBV);
  allG(:,2:ncovp) = gcovBV;
  % combining it with the GP parameters 
  g = - eG' * allG;
  % THE gradient using the EP parameters
  if isFull;
    gcovX = ogpcovgrad(ep.X,gcovX);
    % DERs. wrt. the bias (from ogpcovarp net.bias):
    allEP(:,1) = exp(net.bias)*ones(nAll,1);
    % DERs. wrt. the kernel parameters
    allEP(:,2:ncovp) = gcovX;
    epG = reshape(epG',[nAll,1]);
    g   = - g + epG'*allEP;
  end;
else			      % ONE BY ONE
  % DERs. wrt. the bias (from ogpcovarp net.bias):
  g(1,1) = exp(net.bias)*sum(eG);
  % DERs. wrt. the kernel parameters - USING FOR
  for iKpar = 2:ncovp;
    allG = ogpcovgrad([],gcovBV,iKpar-1);
    g(1,iKpar) = - eG'*allG;
  end;
  % THE gradient wrt. to the left-over parameters.
  if isFull;
    gX(1,1) = exp(net.bias)*sum(epG);
    % DERs. wrt. the kernel parameters - USING FOR
    for iKpar = 2:ncovp;
      allG = ogpcovgrad(ep.X,gcovX,iKpar-1);
      gX(1,iKpar) = epG'*allG;
    end;
    g = - g + gX;
  end;
  fprintf('*');
end;

if find(net.hyplambda);
  pars = (ogppak - net.hypmean);
  g    = g + net.hyplambda.*pars;
end;
