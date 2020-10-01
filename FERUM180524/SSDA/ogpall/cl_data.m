function [X,Y] = cl_data(N,fname,row);
%CL_DATA Generates two-dimensional dataset for classification.
%
%	Description
%
%	[X,Y] = CL_DATA(N,FNAME) - generates data from a mixture of 3
%	Gaussians and labels them with +/-1 where the first class has two
%	components and the second class is a single Gaussian.
%
%	[X,Y] = CL_DATA(N,FNAME,ROW) - if row is nonzero, then the generated
%	inputs are on a grid with the size N (thus the number of examples is
%	NXN) and the corresponding inputs are -1,1.
%
%	See also
%	OGP, OGPTRAIN, C_CLASS_BIN, DEMOGP_REG
%

%	Copyright (c) Lehel Csato (2001-2004)

  
% constants for building the artificial dataset
Disp = -2.5;		      % the shift of the whole dataset
Magnif = 5;		      % magnification factor
% matrix storing the rotations
Rot(:,:,1) = [1 -1;  1 1]/sqrt(2);
Rot(:,:,2) = [1  1; -1 1]/sqrt(2);
Rot(:,:,3) = [1 sqrt(3); -sqrt(3) 1]/2;
L          = [4 2; ...	      % ... the scalings along the diag.s
              4 2; ...
              6 3]/40;
C          = [0.25 0.7; ...   % ... and the individual centres
              0.75 0.7; ...
              0.59 0.35];
pr         = [0.3 0.3 0.4]; % weighting factors
classInd   = [1    1     -1]; % and coefficients

if nargin<3;
  row=0;
end;
if nargin<2;
  fname=[];
end;
    
if row
  X = 0:(1/(N-1)):1;
  [X1, X2] = meshgrid(X);
  X = [X1(:) X2(:)];
  % Computing the class-conditional density
  dM = zeros(size(X,1),2);    % Mahalanobis dist.
  for ii=1:length(pr);
    sSinv = inv(Rot(:,:,ii)*diag(L(ii).^2)*Rot(:,:,ii)');
    cX    = X - repmat(C(ii,:)',1,size(X,1))';
    dM(:,(1.5+classInd(ii)/2)) = dM(:,(1.5+classInd(ii)/2)) + ...
	exp( -1/2*sum((cX*sSinv).*cX,2))*pr(ii);
  end;
  Y  = dM(:,1)./(dM(:,1)+dM(:,2));
else
  % Random mixture appart.
  clRand = rand(N,1);
  clRand = (clRand*ones(size(pr)))> (ones(size(clRand))*cumsum(fliplr(pr)));
  clRand = 3 - sum(clRand,2);
  % Generating inputs
  X      = randn(N,2);
  for ii=1:length(pr);
    iInd = find(clRand==ii);
    X(iInd,:) = ( repmat(C(ii,:)',1,length(iInd))+ ...
		  Rot(:,:,ii)*diag(L(ii,:))*X(iInd,:)')';
  end;
  Y      = classInd(clRand)';
end;

X = Magnif*X + Disp;

if length(fname)
  eval(['save ' fname ' X Y;']);
end;

