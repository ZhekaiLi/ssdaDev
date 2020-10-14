clear; close all; clc
%% demo
EX2;    % call model
global net gpopt ep;


%% Initial function

ferum_ZK

%% case 25: ssda + store GP model in each run to find the prediction with the least variance
% Process
%   For each step,  store the renewed GP model in struct NETS
%   then in prediction process, use all stored GP models to predict
%   and choose the result with the least variance
% Change
%   1. Add a new file called ss_da_compareGPMs.m
%   2. Add a new file called ssda_step_compareGPMs.m


%% case 26: ssda + use Matlab RGP to replace the initial ogp method
% Process
%   Use Matlab's GP to replace the inital ogp (online gaussian process) 
% Change
%   1. Add a new file called ss_da_GPinMatlab.m
%   2. Add a new file called ssda_step_GPinMatlab.m

%% case 27: ssda modified
% Process
%   The case 24 (inital case) use U, G to train but X to predict,
%   this case 27 use X, G to train and X to predict
% Change
%   1. Add a new file called ss_da_modified.m
%   2. Add a new file called ssda_step_modified.m

