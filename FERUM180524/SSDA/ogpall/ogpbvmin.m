function [iMin, scores] = ogpbvmin(gInd);
%OGPBVMIN Finds the BV that contributes the least to the GP 
%
%	Description
%
%	IMIN = OGPBVMIN - returns the index if the BV that has the smallest
%	score (i.e. contribution) to the Gaussian Process NET.  The loss
%	measure is the KL-distance between the GP with the full GP set and
%	the GP with the BV with index IMIN optimally removed.
%
%	If there is a second output variable SCORES, then the loss vector for
%	each BV element is returned.
%
%	If an additional indicator input variable GIND is present, then the
%	KL-distance is replaced with the ``Euclidean'' distance in the
%	feature space, i.e. the values of NET.C and NET.W are ignored.
%
%	See also
%	OGP, OGPTRAIN
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% checking if there is only the geometric part to return
if nargin == 0;
  gInd = 0;
end;

% if there are no BVs, nothing to do
nBV = size(net.BV,1);
if ~nBV
  iMin   = zeros(0,1);
  scores = iMin;
  return;
end;

% if we have one-dimensional outputs, LIFE is simple ...
if net.nout == 1;
  if gInd;
    scores = 1./diag(net.KBinv);
  else
    scores = (net.w.^2)./(diag(net.C) + diag(net.KBinv));
  end;
else			      % for multidim case a lot more work ...
  if gInd;		      % only the inverse determinants
    for iBV = 1:nBV;
      iR = [(1+(iBV-1)*net.nout):(iBV*net.nout)];
      scores(iBV,1) = 1./det(net.KBinv(iR,iR));
    end;
  else
    for iBV = 1:nBV;
      iR = [(1+(iBV-1)*net.nout):(iBV*net.nout)];
      wR = net.w(iR,:);
      scores(iBV,1) = (wR'/(net.KBinv(iR,iR)+net.C(iR,iR)))*wR;
    end;
  end;
end;

% the index of the least 'scorer'
[vMin, iMin] = min(scores);

% BELOW IS AN ALTERNATIVE TO THE COMPUTATION OF SCORES - FOR
% MULTIDIMENSIONAL OUTPUT AND USING NO FOR LOOPS. THE FOR-LOOP SOLUTION
% PREFERRED FOR CLARITY!!

% $$$     indC = [(1-net.nout):(net.nout-1)];
% $$$     mT = spdiags(net.C + net.KBinv,indC);
% $$$     mI = net.nout - 1 + diag(1:size(mT,1))*ones(size(mT)) - ...
% $$$ 	 ones(size(mT))*diag(indC-ceil((indC+1)/2));
% $$$     mI = mod(mI,net.nout) - ones(size(mT))*diag(abs(indC));
% $$$     mT(find(mI<0)) = 0;
% $$$     mT = spdiags(mT,indC,nBV*net.nout,nBV*net.nout);
% $$$     % mT^-1 is copmputed using right division
% $$$     scores = reshape((net.w'/mT).*net.w',[net.nout nBV]);
% $$$     scores = sum(scores,1);
