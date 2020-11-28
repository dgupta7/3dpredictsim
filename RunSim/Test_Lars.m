clear all
close all
clc
%% Paths
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/OCP']);
addpath([pathRepo '/VariousFunctions']);
addpath([pathRepo '/Plots']);

%% Default settings
% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 10;        % number of threads for parallel computing
S.max_iter  = 10000;
solve = 1;
pp = 0;
plot = 0;
% tarsometatarsal joint
S.tmt = 1;
S.tmt_locked = 0;

% exo
S.ExoBool       = 0;
S.ExoScale      = 0;

% output folder
S.ResultsFolder     = 'PredSim_adaptations'; % 'Test_Lars' 'debug_tmt'
S.savename      = 'Pog_s1_tmt_locked_exo_no_ig21'; %'Pog_s1_tmt_d05_k800_no';

% Folder with default functions
S.subject            = 's1_Poggensee';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 1;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)
S.savename_ig   = 'NoExo';



% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% select folder with polynomials
S.PolyFolder = 's1_Poggensee';

% external function
if S.tmt == 0
%     S.CasadiFunc_Folders = 'Casadi_s1Pog_MuscModel_bCst';
    S.CasadiFunc_Folders = 'Casadi_s1Pog_mtp_Default';
    if S.ExoBool*S.ExoScale == 0
        S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtp.dll';        % external function
        S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtp_pp.dll';     % external function for post-processing
    else
        S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp
    end
    
elseif S.tmt ==1
%     S.CasadiFunc_Folders = 'Casadi_s1Pog_tmt_d05_k800';
    S.CasadiFunc_Folders = 'Casadi_s1Pog_tmt_Default';
    if S.ExoBool*S.ExoScale == 0
        S.ExternalFunc  = 'PredSim_3D_Pog_s1_tmt.dll';        % external function
        S.ExternalFunc2 = 'PredSim_3D_Pog_s1_tmt_pp.dll';     % external function for post-processing
    else
        
    end
    
end

% dataset with exoskeleton torque profile
S.DataSet       = 'PoggenSee2020_AFO';





%% Call simulation

if S.tmt == 0
    if solve
    f_PredSim_Gait92(S);        % run the optimization
    end
    if pp
    f_LoadSim_Gait92(S.ResultsFolder,S.savename); % post-proces simulation results
    end
elseif S.tmt == 1
    if solve
    f_PredSim_Gait92_tmt(S);     % run the optimization
    end
    if pp
    f_LoadSim_Gait92_tmt(S.ResultsFolder,S.savename); % post-proces simulation results
    end
end

if plot
    pathResult = fullfile([pathRepo '/Results/' S.ResultsFolder]);
    Plot3D_pwd(pathResult);
end

