% demo
clear; close all; clc;

global net gpopt ep;

EX2;    % call model
star;
gpini(analysisopt,probdata,gfundata,femodel,randomfield);
[ssda_result, probdata ] = ss_da(1,probdata,analysisopt,gfundata,femodel,randomfield);
    Pf_da = [Pf_da ssda_result.pf];
    cov_da = [cov_da ssda_result.cov_pf(1)];
    
analysisopt.echo_flag = 0;
i = 0;
Pf_da = [];
cov_da = [];
nfun_da = [];
accep = {};
while i < 100
    i = i + 1
    [ssda_result, probdata ] = ss_da(1,probdata,analysisopt,gfundata,femodel,randomfield);
    Pf_da = [Pf_da ssda_result.pf];
    cov_da = [cov_da ssda_result.cov_pf(1)];
    y_da{i} = ssda_result.ssda_Data.y;
    p_da{i} = ssda_result.ssda_Data.p;
    nfun_da = [nfun_da ssda_result.nfun];
    accep{i} = ssda_result.ssda_Data.AccRate;
end

% save('Ex2_3.mat','Pf_da','nfun_da','cov_da','accep','y_da','p_da')

j = 0;
Pf_ss = [];
cov_ss = [];
nfun_ss = [];
% accep = {};
while j < 100
    j = j + 1;
    [ss_result, probdata ] = subset_simulation(1,probdata,analysisopt,gfundata,femodel,randomfield);
    Pf_ss = [Pf_ss ss_result.pf];
    cov_ss = [cov_ss ss_result.cov_pf];
    nfun_ss = [nfun_ss ss_result.nfun];
    y_ss{j} = ss_result.SubsetData.y;
    p_ss{j} = ss_result.SubsetData.p;
%     accep{i} = ssda_result.ssda_Data.AccRate;
end

% save('Ex2s_3.mat','Pf_ss','nfun_ss','cov_ss','y_ss','p_ss')

std(Pf_da)/mean(Pf_da)*mean(nfun_da-300)
std(Pf_ss)/mean(Pf_ss)*mean(nfun_ss)
