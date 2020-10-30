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


%% case 260: ssda + use Matlab RGP to replace the initial ogp method
% Process
%   Use Matlab's GP to replace the inital ogp (online gaussian process) 
% Change
%   1. Add a new file called ss_da_GPinMatlab.m
%   2. Add a new file called ssda_step_GPinMatlab.m

%% case 261: case 260 for test (reduce the number of train samples to save time, sacrificing accuracy)
% 以  case 260 为基础，削减了每次高斯过程训练的样本数量（人为使之从1000降到300）
% 进而极大提升训练速度（从60s降到5s）

%% case 262: case 261 modified (modify training samples by finding the nearest base samples with new samples)
% 以 case 261 为基础，对训练集进行了优化
% 使每次生成训练集的方式
% 从取总样本后 300 个样本作为训练集
% 变成从总样本中取出最接近新样本均值的 300 个样本作为训练集

%% case 27: ssda modified
% Process
%   The case 24 (inital case) use U, G to train but X to predict,
%   this case 27 use U, G to train and U to predict
% Change
%   1. Add a new file called ss_da_modified.m
%   2. Add a new file called ssda_step_modified.m



