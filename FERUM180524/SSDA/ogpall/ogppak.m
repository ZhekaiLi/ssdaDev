function hp = ogppak
%
%OGPPAK	Puts the Sparse OGP hyperparametrs into a vector.
%
%	Description
%
%	HP = OGPPAK takes the GLOBAL Gaussian Process data structure NET and
%	combines the hyperparameters into a single row vector HP.  The
%	hyperparameter vector is needed to optimise the OGP structure.
%
%	The transformed hyperparameters returned by this function are used in
%	the MLII optimisation method.
%
%	See also
%	OGP, OGPUNPAK, OGPEVIDGRAD, OGPTRAIN
%

%	Copyright (c) Lehel Csato (2001-2004)


global net gpopt ep;

% Check arguments for consistency
errstring = consist(net, 'ogp');
if ~isempty(errstring);
  error(errstring);
end

hp = [net.bias, net.kpar];
