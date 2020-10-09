clear; close all; clc
%% demo
EX2;    % call model
global net gpopt ep;
global NETS;
NETS = {};


%% Initial function

%ferum

%% ZK: new ferum function that compare the results calculated by different GP models
% For "case 24", 
% the initial ferum.m file use function "ss_da()"
% the modificated ferum_compareGPMs.m file use "ss_da_compareGPMs()"
% ss_da_compareGPMs()  keep the calculated GP model in each step
% In the prediction, the function will do calculation by each GP model
% and choose the results with smallest variation

ferum_compareGPMs



