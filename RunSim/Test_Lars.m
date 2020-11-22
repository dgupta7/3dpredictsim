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
S.v_tgt     = 1.33;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 8;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% Folder with default functions
S.subject            = 's1_Poggensee';

% output folder
S.ResultsFolder     = 'Test_Lars';      % temp folder 

S.CasadiFunc_Folders = 'Casadi_s1Pog_MuscModel_bCst';

% select folder with polynomials
S.PolyFolder = 's1_Poggensee';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';

% external function
% S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp

S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtp.dll';        % external function
S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtp_pp.dll';     % external function for post-processing

% dataset with exoskeleton torque profile
S.DataSet       = 'PoggenSee2020_AFO';

S.ExoBool       = 0;
S.ExoScale      = 0;

S.savename      = 'Test1_bCst_no_v133';

%% Call simulation
% f_PredSim_Gait92(S);     % run the optimization
% f_LoadSim_Gait92(S.ResultsFolder,S.savename); % post-proces simulation results





%% path to scripts, for easy access

open([pathRepo '\CasADiFunctions\CasADiFunctions_Lars.m']);
open([pathRepo '\RunSim\PostProcess_SimulationFolder.m']);
open([pathRepo '\Plots\PlotResultsComparison_3DSim.m']);
open([pathRepo '\Plots\BatchPlotFigs.m']);
open([pathRepo '\VariousFunctions\ModelValidation_CrossCorrelationCoefficient.m']);
open([pathRepo '\VariousFunctions\ModelValidation_NLSError.m']);
open([pathRepo '\RunSim\Validate_Model.m']);
open([pathRepo '\VariousFunctions\footmodel.m']);
open([pathRepo '\VariousFunctions\solveFootmodelParameters.m']);
open([pathRepo '\OCP\f_PredSim_Gait92_tmt.m']);
open([pathRepo '\OCP\f_LoadSim_Gait92_tmt.m']);
open([pathRepo '\Bounds\getBounds_all_tmt.m']);
% open([pathRepo '']);
% open([pathRepo '']);
% open([pathRepo '']);
% open([pathRepo '']);

