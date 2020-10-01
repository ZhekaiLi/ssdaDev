function ogpunpak(hp)
%OGPUNPAK Puts hyperparameters back into the Sparse OGP. 
%
%	Description
%
%	OGPUNPAK(HP) takes the GLOBAL Gaussian process structure NET, the
%	GLOBAL structure EP of TAP/EP parameters together with a
%	hyperparameter vector HP, and returns a Gaussian Process identical to
%	the input model, except that the covariance bias BIAS, the input
%	weight vector parametrised by NET.KPAR(1:NET.NIN) and the covariance
%	function specific parameters in KPAR have all been set to the
%	corresponding elements of HP.  The function updates the inverse
%	kernel matrix.
%
%	The TAP/EP parameters are computed only if we ask for it -- no
%	computation time is wasted.
%
%	See also
%	OGP, OGPPAK, OGPEVIDGRAD, OGPTRAIN
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% Check arguments for consistency
errstring = consist(net, 'ogp');
if ~isempty(errstring); error(errstring); end;

% perform the operations only if the values are different
hp2 = ogppak;
if sum(sum(abs(hp2-hp)))<1e-16; return; end;

wLim = 12;
% log-BIAS of the kernel  function
net.bias = hp(1);
% INPUT scaling
tW = hp(2:(net.nin+1));
tW(find(tW>wLim)) =  wLim;    % approx. 13000
tW(find(tW<-wLim))= -wLim;    % approx. 0.001
% Unpack AMPLITUDE
tA = hp(net.nin+2);
tA(find(tA>wLim))  =  wLim;
tA(find(tA<-wLim)) = -wLim;
% unpack ORDER
tO = hp(net.nin+3);
% special treatment for the POLY case (powers can only be integers).
if strcmp(upper(net.covarfn),'POLY');
  tO = max(log(round(exp(tO))),0);
  tO = min(tO,log(16));
end;

net.kpar = [tW tA tO hp((net.nin+4):end)];

% we need to update the parameters that depend on it:
nBV     = size(net.BV,1);
if gpopt.postopt.isep;
  KBold = net.KB;
end;
net.KB  = ogpcovarp(net.BV,net.BV);
eyeBV   = eye(nBV);
% updating the inverse Gram Matrix via Choleski decomposition
warning off;
net.addEps = 1;
while 1;
  facM          = exp(net.kpar(2)-50+log(net.addEps));
  [net.KBinv,p] = chol(net.KB + eyeBV*facM);
  if ~p | log2(net.addEps)>300;
    break;
  end;
  net.addEps = net.addEps*4;
  fprintf('!');		      % sometimes better WITHOUT warnings
end;
net.KBinv = pinv(net.KBinv);
net.KBinv = net.KBinv*net.KBinv';
% TAP/EP step to change the parameters of
if ~isempty(ep);
  facM     = exp(net.kpar(2)-50+log(net.addEps));
  KB       = net.KB + eyeBV*facM;
  ep.Kplus = ogpcovarp(ep.X,net.BV);
  ep.projP = ep.Kplus/KB;
  UU       = ep.projP' * ep.lamP * ep.projP;
  CC       = UU*KB + eyeBV;
  MM       =   CC \ (ep.projP' * ep.lamP * ep.aP);
  CC       = - CC \ UU;
  net.C = (CC + CC')./2;
  net.w = MM;
else
  matFac = (eyeBV + net.C*(KBold - net.KB));
  net.C  = matFac\net.C;
  net.w  = matFac\net.w;
end;
warning on;
