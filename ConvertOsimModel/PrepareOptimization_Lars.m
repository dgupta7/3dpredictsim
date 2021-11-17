clear
clc
%% Function: prepareInfoNLP

%% Inputs:

% path information
MainPath = 'D:\school\WTK\thesis\model\3dpredictsim';
addpath(fullfile(MainPath,'\VariousFunctions'));
addpath(fullfile(MainPath,'\CasADiFunctions'));
addpath(fullfile(MainPath,'\Polynomials'));
addpath(fullfile(MainPath,'\MetabolicEnergy'));
AddCasadiPaths();



% path to Cpp file used in the optimization
% S.CppFile_NLP = fullfile(MainPath,'ExternalFunctions','CppFiles','PredSim_3D_Fal_s1_v7_test7.cpp');
S.CppFile_NLP = fullfile(MainPath,'ExternalFunctions','CppFiles','SimExo_3D_Pog_s1_mtj_TTC_pp_v1.cpp');

% Number of input arguments in the cpp file
S.CppFile_nInput = 33*3 +2;
% S.CppFile_nInput = 31*3;

% Folder to save the polynomials
S.PolyFolder = 'subject1';
% Modelpath
S.ModelPath = fullfile(MainPath,'OpenSimModel/subject1','subject1_mtj.osim');
% model selection options: Rajagopal, Gait92
S.ModelName = 'Gait92_mtj';   

% path information for automatically building cpp files
OsimSource  = 'D:\opensim-ad\opensim-ad-core';
OsimBuild   = 'D:\opensim-ad\opensim-ad-core-build';
DllPath     = 'D:\school\WTK\thesis\model\3dpredictsim\ExternalFunctions';
ExtFuncs    = 'D:\opensim-ad\external-functions';
VSinstall   = 'C:\Program Files (x86)\Microsoft Visual Studio 14.0';

%% 1: Polynomial fitting

% Fit polynmial functions
% Bool_RunMA = 1; % Boolean to select if we have to run the muscle analysis
% FitPolynomials(MainPath,S.ModelName,S.ModelPath,S.PolyFolder,Bool_RunMA)

%% Create .dll files
% install the functions to create .dll files. you can download this matlab
% software here: https://github.com/MaartenAfschrift/CreateDll_PredSim

% create the .dll file automatically to solve the NLP
[CppDir,Name,~] = fileparts(S.CppFile_NLP);
CreateDllFileFromCpp(CppDir,Name,OsimSource,OsimBuild,DllPath,ExtFuncs,VSinstall,S.CppFile_nInput);

