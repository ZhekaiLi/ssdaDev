function ogpstep_ep(is_sp,indX,K1,K2,gamma,hatE,tA,tV,logZ);
%OGPSTEP_EP Updates the TAP/EP parameters after an online sweep
%
%	Description
%
%	OGPSTEP_EP(IS_SP,INDX,K1,K2,GAMMA,HATE,TA,TV,LOGZ) - updates the
%	GLOBAL structure EP, the parameters of the TAP/EP approximation to
%	reflect the addition of the current input.
%
%	This function is called from OGPPOST and is not used on itw own.
%
%	Depending on the indicator IS_SP, the updates in the structure EP
%	correspond to a full/sparse OGP update.
%
%	See the online (HTML) documentation for details
%
%	See also
%	OGP, OGPTRAIN, OGPSTEP_FULL, OGPSTEP_SP
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% we assume that the update of NET had already been done
iEP = [1+(indX-1)*net.nout:indX*net.nout];

ep.logZ(indX)=logZ+(net.nout*log(2*pi)-log(abs(det(-K2)))-(K1'/K2)*K1)./2;

if is_sp;		      % SPARSE update
  ep.projP(iEP,:)  = hatE';
else			      % FULL update
  nBV          = size(net.BV,1);
  indS         = [1+(nBV-1)*net.nout:net.nout*nBV];
  hatE         = zeros(nBV*net.nout,net.nout);
  hatE(indS,:) = eye(net.nout);
  ep.projP(:,indS)  = 0;
  ep.projP(iEP,:)   = hatE';
end;

ep.aP(iEP)       = tA - K2\K1;
newL = - (eye(net.nout) + K2*tV)\K2;
ep.lamP(iEP,iEP) = (newL + newL')./2;

if is_sp & ~net.proj;	      % if the moment-pres. projection is used
  ep.lamP(iEP,iEP) = ep.lamP(iEP,iEP)/...
      (eye(size(gamma)) + gamma*ep.lamP(iEP,iEP));
end;
