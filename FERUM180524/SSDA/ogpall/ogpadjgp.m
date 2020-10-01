%OGPADJGP Computes the GP coefficients from the TAP/EP ones.
%
%	Description
%
%	OGPADJGP updates the OGP coefficients (W and C) from the GLOBAL
%	TAP/EP structure EP calculated within the online learning routine.
%
%	This function serves as ``sanity check'' to compare the stability of
%	the iterative solution with the ``batch result'' which requires a
%	matrix inversion.  The sum of absolute differences between the new
%	(TAP/EP based) and old (W,C) OGP parameters is displayed.
%
%	This function does not recompute the inverse of the Gram matrix.
%	Recomputing the inverse and adjusting the GP parameters proved to be
%	!EXTREMELY! unstable when DET(NET.KB) is small - which the case in
%	most practical situations when using kernel methods.
%
%	See also
%	OGP, OGPTRAIN
%

%	Copyright (c) Lehel Csato (2001-2004)


global net gpopt ep;

% Recomputing the Gram matrix is STRONLGY unadvised! Probable cause:
%
%     1. the TAP/EP parameters use the iterative approx. to the
%     inverse 
%  AND
%     2. the MATLAB procedure 'inv' or 'pinv' will most probably use a
%     different approximation to the inverse

if ~(isfield(ep,'projP') & isfield(ep,'lamP') & isfield(ep,'aP'));
  error('No valid structure EP');
end;

UU  = ep.projP' * ep.lamP * ep.projP;
CC  = UU*net.KB + eye(size(net.KB));
MM  =   CC \ (ep.projP' * ep.lamP * ep.aP);
CC  = - CC \ UU;

% displaying transformation error
fprintf(['Mean error: ' num2str(sum(abs(net.w-MM)))]);
fprintf([', Cov. error: ' num2str(sum(sum(abs(CC-net.C)))) '\n']);

net.C = CC;
net.w = MM;
