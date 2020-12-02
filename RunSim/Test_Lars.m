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
pp = 0;             % postproces
plot_folder = 0;    % plot figures (for entire resultsfolder)
plot_comp = 0;      % plot figures (for comparison/validation)
timing = 0;         % compare timing of tmt and mtp movement and torque


% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 8;        % number of threads for parallel computing
S.max_iter  = 10;    % maximum number of iterations

% tarsometatarsal joint
S.tmt = 1;              % 1: use a model with tmt joint
S.tmt_locked = 0;       % 1: lock the tmt joint (to compare with model w/o)
S.kTMT = 200;           % (Nm/rad) stiffness of tmt joint 
S.dTMT = 0;           % (Nms/rad) damping of tmt joint

% assumption to simplify Hill-type muscle model
S.MuscModelAsmp = 0;    % 0: musc width = cst, 1: pennation angle = cst

% exo
S.ExoBool       = 1;    % 1: is wearing exo
S.ExoScale      = 0;    % scale factor of exoskeleton assistance profile 
                        % 0: no assistance (passive) 1: nominal assistance (active)
                        
S.DataSet = 'PoggenSee2020_AFO';            % dataset with exoskeleton torque profile
                        
% output folder
S.ResultsFolder = 'PredSim_adaptations';    % other options: 'Test_Lars' 'debug_tmt' 'PredSim_adaptations'


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
pathResults = fullfile([pathRepo '/Results'],S.ResultsFolder);
if ~isfolder(pathResults)
    mkdir(pathResults);
end

% build savename and casadifunction foldername
if strcmp(S.subject,'s1_Poggensee')
    savename = 'Pog_s1';
    casfuncfol = 'casadi_s1Pog';
end
if S.tmt
    savename = [savename '_tmt'];
    casfuncfol = [casfuncfol '_tmt'];
    if S.tmt_locked
        savename = [savename 'L'];
    end
end
if isfield(S,'MuscModelAsmp') && ~isempty(S.MuscModelAsmp) && S.MuscModelAsmp==0
    savename = [savename '_bCst'];
    casfuncfol = [casfuncfol '_MuscModel_bCst'];
else
    savename = [savename '_aCst'];
    casfuncfol = [casfuncfol '_MuscModel_alphaCst'];
end
if S.tmt && isfield(S,'dTMT') && ~isempty(S.dTMT) && ~S.tmt_locked
    savename = [savename '_d0' num2str(S.dTMT*10)];
    casfuncfol = [casfuncfol '_d0' num2str(S.dTMT*10)];
end
if S.tmt && isfield(S,'kTMT') && ~isempty(S.kTMT) && ~S.tmt_locked
    savename = [savename '_k' num2str(S.kTMT)];
    casfuncfol = [casfuncfol '_k' num2str(S.kTMT)];
end
if S.IGsel == 1
    savename = [savename '_ig1'];
else
    savename = [savename '_ig2' num2str(S.IGmodeID )];
end
if S.ExoBool == 1
    if S.ExoScale == 0
        savename = [savename '_pas'];
    else
        savename = [savename '_act'];
    end    
end

S.CasadiFunc_Folders = casfuncfol;
S.savename = savename;

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
        S.ExternalFunc  = 'PredSim_3D_Pog_s1_tmt.dll';        % external function
        S.ExternalFunc2 = 'PredSim_3D_Pog_s1_tmt_pp.dll';     % external function for post-processing
    else
        S.ExternalFunc  = 'SimExo_3D_Pog_s1_tmt.dll';
        S.ExternalFunc2  = 'SimExo_3D_Pog_s1_tmt_pp.dll';
    end
    
end

% Create the casadifunctions if they do not exist yet
if ~isfolder([pathRepo '\CasADiFunctions\' S.CasadiFunc_Folders]) && (solve || pp)
    CreateCasADiFunctions_all_tmt(pathRepo,S);
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

if plot_folder             % plot default figures for entire resultsfolder
    pathResult = fullfile([pathRepo '/Results/' S.ResultsFolder]);
    Plot3D_pwd(pathResult); %'no_meas_data'
end

%%
if plot_comp
    ResultsFile1 = fullfile(pathResults,'Pog_s1_tmt_bCst_d02_k800_ig24_pp.mat');
    ResultsFile2 = fullfile(pathResults,'Pog_s1_tmt_bCst_d02_k1000_ig24_pp.mat');

    
    pathData = [pathRepo,'/ExperimentalData','/ExperimentalData.mat'];
    addpath([pathRepo,'/Plots']);

    % plot default figures
    hh1 = figure();
    set(hh1,'Position',[82         151        1497         827]);
    PlotResults_3DSim_tmt(ResultsFile1,[1 0 0],'k800',hh1);
    PlotResults_3DSim_tmt(ResultsFile2,[0 0 1],'k1000',hh1);
    

    % Compare kinematics and kinetics any number of results to the measurementdata
    ValidationPlots(pathData,ResultsFile1,ResultsFile2);

    % Compare muscles for 2 results
%     PlotResultsComparison_3DSim(ResultsFile1,ResultsFile2);

end

%%

if timing
    ResultsFile = fullfile(pathResults,'Pog_s1_tmt_bCst_d02_k800_ig24_pp.mat');
    Plot_tmt_mtp(ResultsFile)
end

