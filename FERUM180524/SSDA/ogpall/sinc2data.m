function [X,Y] = sinc2data(type,N,var,fname,row,varargin);
%SINC2DATA Generates two-dimensional sinc test data.
%
%	Description
%
%	[X,Y] = SINC2DATA(TYPE,N,VAR,FNAME,ROW) - noisy realisation of
%	Y=A*SINC(|X/D|) where A=10 and the inputs are divided with D=20.
%
%	The parameters:
%
%
%	 TYPE  - the type of the additive noise. Recognised types are:
%	'gauss', 'laplace', and 'posexp'.
%
%	 N  - the number of samples.
%
%	 VAR  - the noise parameter (variance).
%
%	 FNAME  - the name of the file where to write the samples. If the
%	string is empty (or the function has only three parameters), then no
%	file   is written.
%
%	 ROW  - if this option is set, then the input positions X are not
%	random but evenly spaced.
%
%
%	See also
%	OGP, OGPTRAIN, C_REG_GAUSS, C_REG_EXP, C_REG_LAPL, DEMOGP_REG, SINCDATA
%

%	Copyright (c) Lehel Csato (2001-2004)

% Defining default parameters
MagnifY = 5;
MagnifX = 20;

if length(varargin)>0;
  MagnifX = varargin{1};
  if length(varargin)>1; MagnifY = varargin{2}; end
end;

if nargin<4;
  fname= []; row = 0;
elseif nargin<5;
  row  = 0;
end;
if nargin<3;
  error('Not enough arguments in calling SINC2DATA');
end;

  
if row
  nGr      = floor(sqrt(N));
  X        = -1:(2/(nGr-1)):1;
  [xt, yt] = meshgrid(X,X);
  X = [xt(:) yt(:)];
else
  X = 2*rand(N,2)-1;
end;
Y = MagnifY * sinc(sum(X.^2,2));
X = MagnifX * X;

switch upper(type);
 case 'GAUSS'		      % gaussian noise
  Y = Y + sqrt(var)*randn(size(Y));
 case 'LAPLACE'		      % symmetric laplace noise
  if var
    Y = Y + laplace(size(Y))/var;
  end;
 case 'POSEXP'		      % nonsymmetric exp. noise
  if var;
    Y = Y + laplace([N,1],0)./var;
  end;
 otherwise
  error('Unknown noise type');
end;

if length(fname)
  eval(['save ' fname ' X Y;']);
end;
