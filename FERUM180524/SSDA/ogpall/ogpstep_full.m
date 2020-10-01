function ogpstep_full(cX,kX,kk,K1,K2,gamma,hatE);
%OGPSTEP_FULL Performs a full online update step of the Gaussian Process
%
%	Description
%
%	[NET] = OGPSTEP_FULL(NET,CX,INDX,KX,KK,K1,K2,GAMMA,HATE) - updates
%	the GLOBAL Gaussian process structure NET with the information from
%	the current data point.  This function is called from OGPPOST and
%	this function makes no change inthe EP-representation of the GP (the
%	chenges are done in OGPSTEP_EP).         The BV set is extended with
%	CX - the current input example. GAMMA and HATE are the geometric
%	components required to update the inverse of the Gram, or kernel
%	matrix.
%
%	See the online (HTML) documentation for details
%
%	See also
%	OGP, OGPTRAIN, OGPSTEP_SP, OGPSTEP_EP
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% computing index ranges
nBV  = size(net.BV,1) + 1;
indS = [(1+(nBV-1)*net.nout):(nBV*net.nout)];
indL = [ 1:((nBV-1)*net.nout)];

% updating KBinv
net.BV(nBV,:)        = cX;
hatE(indS,:)         = -eye(net.nout);
net.KBinv(indS,indS) = 0; 
net.KBinv            = net.KBinv + hatE*inv(gamma)*hatE';

% updating the Gram matrix
net.KB(indS,indS) = 0;
net.KB(indL,indS) = kX;
net.KB(indS,indL) = kX';
net.KB(indS,indS) = kk;

% updating GP parameters
tVector          = net.C*kX;
tVector(indS,:)  = eye(net.nout);
% extending the size of the parameters
net.w(indS,:)    = 0; 
net.C(indS,indS) = 0;

net.w = net.w + tVector*K1;
net.C = net.C + tVector*K2*tVector';

% stability
net.C     = (net.C  + net.C' )./2;
net.KB    = (net.KB + net.KB')./2;
net.KBinv = (net.KBinv + net.KBinv')./2;
