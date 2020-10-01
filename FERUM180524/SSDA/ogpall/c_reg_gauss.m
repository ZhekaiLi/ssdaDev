function [logLik,K1,K2,varargout] = c_reg_gauss(likpar,y,mX,sigX2,varargin)
%C_REG_GAUSS Online update coeffs. for regression with Gaussian noise.
%
%	Description
%
%	[LOGLIK,K1,K2] = C_REG_GAUSS(LIKPAR,Y,MX,VARX) - takes the noise
%	variance (likelihood parameter) together with the mean, variance and
%	the prior mean at the current input and returns the coefficients for
%	the online update of the GP quadratic regression.
%
%	Parameters:
%
%	 LIKPAR  - the noise variance (from NET.LIKPAR).
%
%	 Y  - observed noisy output.
%
%	 MX  - mean of the GP marginal at the new input.
%
%	 VARX  - the variance of the GP marginal at the new input.
%
%	 LOGLIK  - the logarithm of the averaged likelihood function.
%
%	 K1,K2  - the first and second derivatives of LOGLIK.
%
%	If there are more input points, then all outputs K1,K2,LOGLIK will be
%	matrices [N,1] where N is the numober of observations (columns) in
%	the inputs.
%
%	See also
%	OGP, OGPTRAIN, ERR_MSE, C_REG_LAPL, C_CLASS_BIN, DEMOGP_REG
%

%	Copyright (c) Lehel Csato (2001-2004)

varargout = {[]};

nout = size(mX,2);
nX   = size(y,1);
nM   = size(mX,1);

if size(mX,1) ~= size(sigX2,1)./nout;
  error('Lengths of mean and (co)variance terms do not match\n\n');
elseif nX>1 && nM>1 && nX~=nM;
  error('Lengths of means and observations should match');
end;

if length(likpar)==1;
  covN = likpar*eye(nout);
else
  covN = diag(likpar);
end;

% depending on the dimension of output ...
if nout>1;
  sig2C = sigX2;
  mC    = mX;
  yC    = y;
  for iX = 1:nX;	      % doing sequentially
    if nM>1;
      sig2C = sigX2(1:nout,(nout*(iX-1)+1):(nout*iX));
      mC    = mX((nout*(iX-1)+1):(nout*iX),1);
    end;
    if nX>1;
      yC = y((nout*(iX-1)+1):(nout*iX),1);
    end;
    % checking the validity of the variance for the GP marginal
    if det(sig2C) < 0;
      warning(['Pred. distr. var.: ' num2str(det(sigX2))]);
      sig2C = sig2C - (det(sigX2)+0.05)*eye(nout);
    end;
    % derivative of the log-average
    sig2C      = sig2C + covN;
    K2(iX)     = - sig2C\eye(nout);
    K1(iX)     = - K2*(yC - mC);
    logLik(iX) = - ( log(det(2*pi*sig2C)) + (yC-mC)'*K1)./2;
  end;
else			      % one-dimensional output
  sigX2  = sigX2 + covN;
  K2     = - 1./sigX2;
  K1     = - K2 .* (y-mX);
  logLik = - (log(2*pi*sigX2) + (y-mX).*K1)./2;
end;
