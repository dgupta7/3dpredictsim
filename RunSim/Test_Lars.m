clear all
close all
clc

%% Paths
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/OCP']);
addpath([pathRepo '/VariousFunctions']);
addpath([pathRepo '/Plots']);
addpath([pathRepo '/CasADiFunctions']);
addpath([pathRepo '/Musclemodel']);
addpath([pathRepo '/Polynomials']);

%% Manual settings
solve = 1;          % run solver
pp = 1;             % postproces
plot = 0;      % plot solution

% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 8;        % number of threads for parallel computing
S.max_iter  = 10000;    % maximum number of iterations

% tarsometatarsal joint
S.tmt = 0;              % 1: use a model with tmt joint
S.tmt_locked = 0;       % 1: lock the tmt joint (to compare with model w/o)
S.kTMT = 2000;           % (Nm/rad) stiffness of tmt joint 
S.dTMT = 0;           % (Nms/rad) damping of tmt joint
% nonlinear spring tmt
S.TMT_linear = 1;
S.k1TMT = 800;
S.k2TMT = 1;
S.t1TMT = 0.5;

% assumption to simplify Hill-type muscle model
S.MuscModelAsmp = 0;    % 0: musc width = cst, 1: pennation angle = cst

% exo
S.ExoBool       = 0;    % 1: is wearing exo
S.ExoScale      = 0;    % scale factor of exoskeleton assistance profile 
                        % 0: no assistance (passive) 1: nominal assistance (active)
                        
S.DataSet = 'PoggenSee2020_AFO';            % dataset with exoskeleton torque profile
                        
% output folder
S.ResultsFolder = 'debug';


% Folder with default functions
S.subject            = 's1_Poggensee';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)

if S.IGmodeID == 4
S.savename_ig   = 'NoExo';
elseif S.IGmodeID == 3
S.ResultsF_ig   = 'PredSim_adaptations';
S.savename_ig   = '';
end


%% Automated settings

% make folder to store results if it doesn't exist
pathResults = fullfile([pathRepo '/Results'],S.ResultsFolder);
if ~isfolder(pathResults)
    mkdir(pathResults);
end

% build standardised names
[savename, casfuncfol] = getSavename(S);
S.CasadiFunc_Folders = [casfuncfol '_test'];
S.savename = [savename '_test'];

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% select folder with polynomials
S.PolyFolder = 's1_Poggensee';

% external function
if S.tmt == 0
    if S.ExoBool == 0
        S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtp.dll';        % external function
        S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtp_pp.dll';     % external function for post-processing
    else
        S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp
    end
    
elseif S.tmt ==1
    if S.ExoBool == 0
        S.ExternalFunc  = 'PredSim_3D_Pog_s1_tmt_v3.dll';        % external function
        S.ExternalFunc2 = 'PredSim_3D_Pog_s1_tmt_pp_v3.dll';     % external function for post-processing
    else
        S.ExternalFunc  = 'SimExo_3D_Pog_s1_tmt_v4.dll';
        S.ExternalFunc2  = 'SimExo_3D_Pog_s1_tmt_pp_v4.dll';
    end
    
end

% Create the casadifunctions if they do not exist yet
if ~isfolder([pathRepo '\CasADiFunctions\' S.CasadiFunc_Folders])
    disp('Creating casadifunctions...');
    CreateCasADiFunctions_all_tmt(pathRepo,S);
    disp('...casadifunctions created');
end


%% Run

if S.tmt == 0
    if solve        % run the optimization
    f_PredSim_Gait92(S);
    end
    if pp           % post-proces simulation results
    f_LoadSim_Gait92(S.ResultsFolder,S.savename);
    end
elseif S.tmt == 1
    if solve        % run the optimization
    f_PredSim_Gait92_tmt(S);
    end
    if pp           % post-proces simulation results
    f_LoadSim_Gait92_tmt(S.ResultsFolder,S.savename);
    end
end



%%
if plot
    h = figure();
    set(h,'Position',[82 151 1497 827]);
    PlotResults_3DSim_tmt(fullfile(pathResults,[S.savename '_pp.mat']),[1 0 0],savename,h);
end

% fullfile(pathResults,[S.savename '_pp.mat'])


