clear 
close all
clc


%% Paths

[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/VariousFunctions']);
dpath = [pathRepo,'/Results'];
pathData = [pathRepo,'/ExperimentalData','/ExperimentalData.mat'];

ResultsFile1 = fullfile(dpath,'Test_Lars','Test1_bCst_no_pp.mat');
ResultsFile2 = fullfile(dpath,'Test_Lars','Test1_alphaCst_no_pp.mat');


% ccc is the cross-correlation coefficient of the simulation data and the
% mean measurement data

[rmse1,ccc1] = ModelValidation(ResultsFile1, pathData);
[rmse2,ccc2] = ModelValidation(ResultsFile2, pathData);

rmse = [rmse1,rmse2]
ccc = [ccc1,ccc2]

