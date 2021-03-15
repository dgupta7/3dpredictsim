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
% S.CppFile_NLP = fullfile(MainPath,'ExternalFunctions','CppFiles','PredSim_3D_Fal_s1_tmt_pp_v1.cpp'); 
% S.CppFile_NLP = fullfile(MainPath,'ExternalFunctions','CppFiles','Foot_3D_Pog_s1_mtj_v4.cpp');
S.CppFile_NLP = fullfile(MainPath,'ExternalFunctions','CppFiles','Foot_3D_Ker_87_mtj_v4.cpp'); 

% Number of input arguments in the cpp file
S.CppFile_nInput = 22*3; 


% path information for automatically building cpp files
OsimSource  = 'D:\opensim-ad\opensim-ad-core';
OsimBuild   = 'D:\opensim-ad\opensim-ad-core-build';
DllPath     = 'D:\school\WTK\thesis\model\3dpredictsim\ExternalFunctions';
ExtFuncs    = 'D:\opensim-ad\external-functions';
VSinstall   = 'C:\Program Files (x86)\Microsoft Visual Studio 14.0';


%% 3: Create .dll files if needed

% install the functions to create .dll files. you can download this matlab
% software here: https://github.com/MaartenAfschrift/CreateDll_PredSim
% addpath('C:\Users\u0088756\Documents\FWO\Software\GitProjects\CreateDll_PredSim');

% create the .dll file automatically to solve the NLP
[CppDir,Name,~] = fileparts(S.CppFile_NLP);
CreateDllFileFromCpp(CppDir,Name,OsimSource,OsimBuild,DllPath,ExtFuncs,VSinstall,S.CppFile_nInput);

