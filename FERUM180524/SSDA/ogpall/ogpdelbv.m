function ogpdelbv(iBV);
%OGPDELBV Deletes the specified BVs from the BV set of the GP.
%
%	Description
%
%	OGPDELBV(IBV) deletes the given BVs from the set of BVs.
%
%	The operations are performed on the elements of the structures NET
%	and EP.
%
%	If the TAP/EP parameters have also been kept, then -- together with
%	the removal of the BV's from the NET structure -- the respective
%	parts of the structure EP are also updated.
%
%	Parameters:
%
%	 IBV  - will be deleted. iBV is a row of indices.
%
%	 NET  - the GLOBAL OGP structure from which the BVs with index
%
%	 EP  - the GLOBAL structure storing the TAP/EP parameters.
%
%	See also
%	OGP, OGPTRAIN
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

%%%%%%%%%%%%%%%% ARGUMENT CHECKING %%%%%%%%%%%%%%%%%%
nBV = size(net.BV,1);
if size(iBV,1) ~= 1;
  error('Wrong format for the BV indices');
elseif max(iBV)>nBV | min(iBV)<1;
  error('Indices out of range');
end;
if nargout==2 & nargin~=3;
  error('EP must appear as the third input argument');
end;

%%%%%%%%%%%%%%%% REMOVING A BASIS VECTOR %%%%%%%%%%%%
% if there is no inverse Gram matrix, REcompute
if size(net.KBinv,2)~=nBV*net.nout;
  net.KBinv = ogpcovarp(net.BV,net.BV);
  net.KBinv = pinv(net.KBinv);
  warning('Inverse covariance recomputed!');
end;

% computing indices
rmI = net.nout*repmat((iBV-1)',[1,net.nout]) + ...
      ones(length(iBV),net.nout)*diag([1:net.nout]);
rmI = rmI'; rmI = rmI(:)';
exI = setdiff([1:(nBV*net.nout)],rmI);
exS = setdiff([1:nBV],iBV);

% Computing auxilliary quantities
q_star = net.KBinv(rmI,rmI);
red_q  = net.KBinv(exI,rmI);
c_star = net.C(rmI,rmI);
red_Ag = red_q + net.C(exI,rmI);

% updating GP params
if net.proj;		      % update that keeps the can.pars.
  net.w  = net.w(exI,:) - red_Ag * ((q_star+c_star)\net.w(rmI));
  net.C  = net.C(exI,exI) + red_q * (q_star\red_q') ...
	   - red_Ag * ((q_star+c_star)\red_Ag');
else			      % update that KEEPS the marginals
  tempQ = red_q/q_star;
  net.w = net.w(exI) - tempQ*net.w(rmI);
  red_c = net.C(rmI,exI);
  net.C = net.C(exI,exI) + tempQ*c_star*tempQ';
  tempQ = tempQ*red_c;
  net.C = net.C - tempQ - tempQ';
end;
% stability
net.C     = (net.C  + net.C')./2;

if ~isempty(ep) & isfield(ep,'lamP'); % updating EP parameters
  % Additoinal operations needed when the projection is the MOMENT-MATCHING
  % ONE!
  if ~net.proj;
    if length(iBV)>1;
      error('Basis Vectors have to be removed one by one');
    end;
    % We have to transform 'IBV' - indices in the BV set into the index in
    % the training set.
    prDel = ep.projP(:,iBV).^2;
    lDel  = full(ep.lamP*prDel);
    ep.lamP = ep.lamP + diag(lDel.^2./(q_star + lDel));
  end;
  ep.projP = ep.projP(:,exI) - ep.projP(:,rmI)*(q_star\red_q');
end;

% updating invGram
net.KBinv  = net.KBinv(exI,exI) - red_q*(q_star\red_q');
net.KBinv  = (net.KBinv + net.KBinv')./2;

% updating the Gram matrix
net.KB   = net.KB(exI,exI);
net.BV   = net.BV(exS,:);
