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

% settings:
% Folder to save the polynomials
S.PolyFolder = 's1_Poggensee';
% Modelpath
S.ModelPath = fullfile(MainPath,'OpenSimModel','Subject1_Poggensee.osim'); 
% Folder with CasadiFunctions
S.CasadiFunc_Folders = 'debug_notmt'; %'Casadi_s1Pog_tmt_d05_k800'; 
% path to Cpp file used in the optimization
S.CppFile_NLP = fullfile(MainPath,'ExternalFunctions','CppFiles','PredSim_3D_Pog_s1_tmt_v2.cpp'); 
 % path to Cpp file for post processing
S.CppFile_pp = fullfile(MainPath,'ExternalFunctions','CppFiles','PredSim_3D_Pog_s1_tmt_pp_v2.cpp');    
% Number of input arguments in the cpp file
S.CppFile_nInput = 33*3; 
% model selection options: Rajagopal, Gait92
S.ModelName = 'Gait92';      

% specific settings for exporting casadi functions
% SettingsCasFunc.kTendon_CalfM = 20;
% SettingsCasFunc.kMTP = 1.5/(pi/180)/5;
% SettingsCasFunc.dMTP = 0.5;

SettingsCasFunc.tmt = 0;
% SettingsCasFunc.kTMT = 800; % 250, 500, 1000,2000, 4000
% SettingsCasFunc.dTMT = 0.5;


% path information for automatically building cpp files
OsimSource  = 'D:\opensim-ad\opensim-ad-core';
OsimBuild   = 'D:\opensim-ad\opensim-ad-core-build';
DllPath     = 'D:\school\WTK\thesis\model\3dpredictsim\ExternalFunctions';
ExtFuncs    = 'D:\opensim-ad\external-functions';
VSinstall   = 'C:\Program Files (x86)\Microsoft Visual Studio 14.0';

%% 1: Polynomial fitting

% % Fit polynmial functions
% Bool_RunMA = 0; % Boolean to select if we have to run the muscle analysis
% FitPolynomials(MainPath,S.ModelName,S.ModelPath,S.PolyFolder,Bool_RunMA)

%% 2: Create Casadifunctions

% create casadi functions for equations in optimiztion problem
% CreateCasadiFunctions(MainPath, S.ModelName, S.ModelPath, S.CasadiFunc_Folders,S.PolyFolder,SettingsCasFunc);
   

%% 3: Create .dll files if needed

% install the functions to create .dll files. you can download this matlab
% software here: https://github.com/MaartenAfschrift/CreateDll_PredSim
% addpath('C:\Users\u0088756\Documents\FWO\Software\GitProjects\CreateDll_PredSim');


% % create the .dll file automatically to solve the NLP
[CppDir,Name,~] = fileparts(S.CppFile_NLP);
CreateDllFileFromCpp(CppDir,Name,OsimSource,OsimBuild,DllPath,ExtFuncs,VSinstall,S.CppFile_nInput);

% % create the .dll file automatically for postprocessing
[CppDir,Name,~] = fileparts(S.CppFile_pp);
CreateDllFileFromCpp(CppDir,Name,OsimSource,OsimBuild,DllPath,ExtFuncs,VSinstall,S.CppFile_nInput);
