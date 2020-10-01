%OGPRESET Resets the Gaussian Process.
%
%	Description
%
%	OGPRESET resets the GLOBAL Gaussian Process structure NET to its
%	original state.
%
%	If elements can be added to the BV set, then the BV set itself is
%	cleared, otherwise (ie. if NET.ISBVFIXED~=0) the BV set is kept but
%	other Gaussian Process parameters in NET are set to zero.  One can
%	use this function to reset the OGP if kernel parameters change.
%
%	If there exists a nonempty GLOBAL structure EP, then the fields in
%	this structure are also reset to a default value.
%
%	See also
%	OGP, OGPINIT, OGPTRAIN
%

%	Copyright (c) Lehel Csato (2001-2004)


global net gpopt ep;

% remove the BV set only if we can add new!
if isfield(ep,'X');
  tDim = size(ep.X,1)*net.nout;
else
  tDim = 0;
end;

if ~(net.isBVfixed)
  net.BV    = zeros(0,net.nin);
  net.KB    = zeros(0,0);
  net.KBinv = zeros(0,0);
  net.w     = zeros(0,1);
  net.C     = zeros(0,0);
else
  bvDim     = size(net.BV,1)*net.nout;
  net.w     = zeros(bvDim,1);
  net.C     = zeros(bvDim,bvDim);
end;
ep = [];
