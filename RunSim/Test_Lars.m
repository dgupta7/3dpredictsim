clear all
close all
clc

%% Paths

[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/OCP']);
addpath([pathRepo '/VariousFunctions']);


%% Default settings

% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 2;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% Folder with default functions
S.subject            = 's1_Poggensee';

% output folder
S.ResultsFolder     = 'Test_Lars';      % temp folder 

S.CasadiFunc_Folders = 'Casadi_s1Pog_ScaleParam_k35_bCst';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';

% external function
S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp

% dataset with exoskeleton torque profile
S.DataSet       = 'PoggenSee2020_AFO';

S.ExoBool       = 1;
S.ExoScale      = 1;

S.savename      = 'Test1_exo';

%% Call simulation
f_PredSim_PoggenSee2020(S);





%% path to scripts, for easy access

% D:\Gebruiker\Documents\~school~\WTK\thesis\model\3dpredictsim_Lars\CasADiFunctions\CasADiFunctions_Lars.m
% D:\Gebruiker\Documents\~school~\WTK\thesis\model\3dpredictsim_Lars\RunSim\PostProcess_SimulationFolder.m


