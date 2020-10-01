%% demo
EX1;    % call model
% global net gpopt ep;

ssini;  % initialize reliability problem
gpini(analysisopt,probdata,gfundata,femodel,randomfield);   % initialize gp model
[ssda_result, probdata ] = ss_da(1,probdata,analysisopt,gfundata,...
    femodel,randomfield);  % main program


