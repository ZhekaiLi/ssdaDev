function ogpemptybv(cX,iOpt,addMax);
%OGPEMPTYBV Adds input elements to the BV set without altering the GP.
%
%	Description
%
%	OGPEMPTYBV(CX) - increases the BV set in the Gaussian Process NET by
%	adding the inputs from CX.
%
%	The addition to the BV set is sequential (i.e. a single element of CX
%	at each step) and if the feature-space distance GAMMA between the
%	subspace determined by the BV set and the input is smaller than the
%	threshold (stored in the structure NET), then the respective input is
%	not added.  This helps mantaining a well-conditioned Gram matrix.
%
%	A second form of adding elements to the BV set is to add them
%	according to the "novelty" they bring, measured by the distance from
%	the subspace determined by the existing BVs. This requires a nonzero
%	value IOPT as the third argument to the function:
%	[net] = ogpemptybv(net,cx,iOpt)
%
%	There is a possibility to limit the number of added BVs by specifying
%	a fourth argument
%
%	[net] = ogpemptybv(net,cx,iOpt,addMax)
%
%	which limits the maximum number of new BV additions to ADDMAX.
%
%	If you want to limit the number of BVs and not have addition based on
%	the distances (ie. IOPT==0), then you should call the function with
%	four arguments, where the third is set to zero.
%
%	The additional output argument RESIDERR tells the maximal distance of
%	the new subspace from the input elements that were not included into
%	the new BV set.
%
%	See also
%	OGP, OGPINIT, OGPTRAIN
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% Check arguments for consistency
errstring = consist(net, 'ogp', cX);
if ~isempty(errstring);
  error(errstring);
end;

% setting default argument iOpt to 0 if it is not defined
if nargin==1;
  iOpt   = 0;
  addMax = size(cX,1);
end;
if nargin==2;
  addMax = size(cX,1);
end;

resid = 0;
nBV   = size(net.BV,1);
maxBV = nBV + addMax;
% checking whether we have a valid inverse Gram matrix
if size(net.KBinv,2) ~= nBV*net.nout;
  net.KB    = ogpcovarp(net.BV,net.BV);
  net.KBinv = net.KB\eye(size(net.KB));
  warning('Warning: inverse covariance recomputed!');
end;

for iNew = 1:size(cX,1);
  % quitting if the maximum number of iterations has been exceeded
  if size(net.BV,1)>=maxBV;
    return;
  end;
  % finding the most relevant element in cX
  if iOpt;
    % computing distances
    kX = ogpcovarp(net.BV,cX);
    kk = ogpcovarp(cX);
    sc = -ones(size(cX,1),1);
    for iX = 1:size(cX,1);
      tInd   = [(1+net.nout*(iX-1)):(net.nout*iX)];
      tK     = kX(:,tInd);
      sc(iX) = det(kk(tInd,:) - tK'*(net.KB\tK));
    end;
    [resid,iM]  = max(sc);
    tVar     = cX(1,:);
    cX(1,:)  = cX(iM,:);
    cX(iM,:) = tVar;
    tInd     = [(1+net.nout*(iM-1)):(net.nout*iM)];
    kX       = kX(:,tInd);
    kk       = kk(tInd,:);
    if resid<net.thresh;      % if small, no more additions
      return;
    end;
  else
    % computing the kernel products
    kX = ogpcovarp(net.BV,cX(1,:));
    kk = ogpcovarp(cX(1,:));
  end;
  % computing the geom. components
  hatE  = net.KB\kX;
  gamma = kk - hatE'*kX;
  % adding the new BV
  if det(gamma)>net.thresh;
    indS    = [1:net.nout] + nBV*net.nout;
    indL    = [1:net.nout*nBV];
    nBV     = nBV + 1;
    % updating KBinv
    net.BV(nBV,:)        = cX(1,:);
    hatE(indS,:)         = - eye(net.nout);
    net.KBinv(indS,indS) = 0; 
    net.KBinv            = net.KBinv + hatE*(gamma\hatE');
    % net.KBinv            = (net.KBinv + net.KBinv')/2;
    net.KB(indS,indS)    = 0;
    net.KB(indL,indS)    = kX;
    net.KB(indS,indL)    = kX';
    net.KB(indS,indS)    = kk;
    % updating GP
    net.w(indS,1)    = 0;
    net.C(indS,indS) = 0;
    % updating residual measurements
    resid = det(gamma);
  elseif resid<gamma;
    resid = det(gamma);
  end;
  % removing the first row of cX - the processed element 
  lX  = size(cX,1);
  cX  = cX(2:lX,:);
end;
