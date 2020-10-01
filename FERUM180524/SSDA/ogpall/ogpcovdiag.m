function covf = ogpcovdiag(x1)
%OGPCOVDIAG Calculates the diagonal of the covariance function for the OGP.
%
%	Description
%
%	[COVF] = OGPCOVDIAG(NET, X1) takes an OGP data structure NET together
%	with a data matrix X1 and computes the diagonal elements of the
%	covariance function.
%
%	If the output is not one-dimensional, then the returned matrix is
%	[DIM(X1),NOUT,NOUT] dimensional.
%
%	See also
%	OGP, OGPINIT, OGPCOVARF, OGPCOVARP
%

%	Copyright (c) Lehel Csato (2001-2004)

global net gpopt ep;

% consistency checking
errstring = consist(net, 'ogp', x1);
if ~isempty(errstring);
  error(errstring);
end

beta = spdiags(exp(net.kpar(1:net.nin))',...
	       0, ...
	       net.nin*net.nout,net.nin*net.nout);
aA   = exp(net.kpar(net.nin+1));
ord  = exp(net.kpar(net.nin+2));

n1   = size(x1,1); 	      % OGPCOVARF(X1,X1)
covf = repmat(eye(net.nout),[n1 1])*0;
switch upper(net.covarfn); 
 case {'SQEXP', 'RATQUAD', 'POLYEXP', 'OU'}
  covf = aA*exp(covf);
 case 'POLY';
  covf = aA * (1 + sum((x1*beta).*x1,2) ).^round(ord);
 case 'LINSPLINE'
  if find(x1(:) < 0);
    error('Nodes for LINSPLINE need to be positive');
  else
    for iIn=1:net.nin;
      covf = covf + beta(iIn,iIn)*(x1(:,iIn).^3/3 + x1(:,iIn).^2 +1);
    end;
  end;
  
 case 'SINC'
  tPar = exp(net.kpar((net.nin+1):end));
  if length(tPar==1);
    tPar(2:net.nin) = tPar;
  end;
  for iIn = 1:net.nin;
    covf = covf + tPar(iIn)*ones(n1,1);
  end;
  
 case 'FSIN';		      % FINITE dimensional periodic kernel
  for iIn = 1:net.nin;
    covf = covf + aA * (2*ord+1) * ones(n1,1);
  end;

 case 'USER';		      % user-specified covariance funtion.
  covf = feval(net.kfnaddr,x1);

 otherwise
  error(['Unknown covariance function ', net.covarfn]);
end;
