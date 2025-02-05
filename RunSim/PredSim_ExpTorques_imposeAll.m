
%% Imposing stride frequency
%-----------------------------

% Influence of the kinematic constraint on stepwidth
clear all; close all; clc;

% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 4;        % number of threads for parallel computing

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% Folder with default functions
S.subject            = 's1_Poggensee';

% output folder
S.ResultsFolder     = 'ExpTorques_SW_Freq_Timing2';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';

% select the CasadiFolder
S.CasadiFunc_Folders = 'Casadi_s1Pog_mtp';

% stride frequency between 0.99 and 1.01
S.Bounds.tf = [0.495 0.505];

% Timing exoskeleton assistance
S.PercStance.bool = 1;
S.PercStance.xStanceOr = 0.61;
S.PercStance.xStanceNew = 0.58;

% normal walking simulation
S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtp.dll';        % this one is with the pinjoint mtp
S.ExternalFunc2  = 'PredSim_3D_Pog_s1_mtp_pp.dll';        % this one is with the pinjoint mtp
S.ExoBool       = 0;    
S.ExoScale      = 0;
S.savename      = 'NoExo';
SWadd = 0.13;
S.Constr.calcn = 0.09 + SWadd;
S.Constr.toes = 0.1 + SWadd;
S.Constr.tibia = 0.11 + SWadd;
f_PredSim_PoggenSee2020(S);

% passive simulation
S.DataSet       = 'PoggenSee2020_ExpPass';
S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp
S.ExternalFunc2 = 'SimExo_3D_ExportAll.dll';        % external function for post-processing
S.ExoBool       = 1;    
S.ExoScale      = 1;
S.savename      = 'Passive';
SWadd = 0.15;
S.Constr.calcn = 0.09 + SWadd;
S.Constr.toes = 0.1 + SWadd;
S.Constr.tibia = 0.11 + SWadd;
f_PredSim_PoggenSee2020(S);

% active simulation
S.DataSet       = 'PoggenSee2020_Exp';
S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp
S.ExternalFunc2 = 'SimExo_3D_ExportAll.dll';        % external function for post-processing
S.ExoBool       = 1;    
S.ExoScale      = 1;
S.savename      = 'Active';
SWadd = 0.15;
S.Constr.calcn = 0.09 + SWadd;
S.Constr.toes = 0.1 + SWadd;
S.Constr.tibia = 0.11 + SWadd;
f_PredSim_PoggenSee2020(S);
