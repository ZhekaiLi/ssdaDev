clear; close all; clc
%% demo
EX2;    % call model
global net gpopt ep;


%% Initial function

%ferum

%% ZK: new ferum function that compare the results calculated by different GP models
% For "case 24", 
% the initial ferum.m file use function "ss_da()"
% the modificated ferum_compareGPMs.m file use "ss_da_compareGPMs()"
% ss_da_compareGPMs()  keep the calculated GP model in each step
% In the prediction, the function will do calculation by each GP model
% and choose the results with smallest variation

% global NETs;
% NETs = {};
%ferum_compareGPMs

%% ZK: new ferum function that use different GP in matlab, and compare the results calculated by different GP models
% For "case 24", 
% the initial ferum.m file use function "ss_da()"
% the modificated ferum_GPinMatlab.m file use "ss_da_GPinMatlab()"
% ss_da_GPinMatlab()  keep the calculated GP model in each step
% In the prediction, the function will do calculation by each GP model
% and choose the results with smallest variation

global data;
global gpm;
ferum_GPinMatlab

