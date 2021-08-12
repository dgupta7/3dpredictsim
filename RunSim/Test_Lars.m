%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script serves as the main file for the changes made in the branch 
% Thesis_Lars. It allows to specify the settings, solve, post-process and
% plot a single gait simulation.
%
% Alternatively, the specified settings can be saved and added to a batch.
% To run the batch, run \BatchRunQueue.m
%
% This script can also be used to run a static simulation of a single foot.
% This foot can be hanging in the air supported at the knee, or standing on
% the ground with a vertical force pushing down on the knee.
%
% To plot multiple results on the same figure, use \Plots\make_any_plot.m
%
% Author: Lars D'Hondt (March 2021)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
% close all
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
addpath([pathRepo '/Debug']);
addpath([pathRepo '/FootModel']);
AddCasadiPaths();

%% General settings
% Full body gait simulation
slv = 1;                % run solver
pp =1;                 % postproces
plot = 1;               % plot solution
batchQueue = 0;         % save settings to run later
% Static foot simulation
foot_standing = 0;      % load on knee
foot_hanging = 0;       % knee position fixed, tibia and foot hanging freely

% settings for optimization
S.v_tgt     = 1.25;     % average speed 1.25
S.N         = 50;       % number of mesh intervals
S.NThreads  = 6;        % number of threads for parallel computing
% S.max_iter  = 10;       % maximum number of iterations (comment -> 10000)

% output folder
S.ResultsFolder = 'MidTarsalJoint'; % 'MuscleModel' 'MidTarsalJoint' 'Final'
suffixCasName = '_PFx10';         
suffixName = '_PFx10';

% assumption to simplify Hill-type muscle model
S.MuscModelAsmp = 0;    % 0: musc width = cst, 1: pennation angle = cst
S.contactStiff = 10;    % 10: contact spheres are 10x stiffer

% Test subject
S.subject = 'subject1'; % barefoot
% S.subject = 's1_Poggensee'; % normal shoes, exoskeleton powered on/off

%% Tarsometatarsal joint
S.tmt = 0;              % 1: use a model with tmt joint
% linear spring
S.kTMT = 1000;          % (Nm/rad) stiffness of tmt joint 
S.dTMT = 0.5;           % (Nms/rad) damping of tmt joint
% nonlinear spring tmt
S.TMT_linear = 1;
S.k1TMT = 800;
S.k2TMT = 1;
S.t1TMT = 0.5;
% windlass mechanism
S.Windlass = 1;
S.cWL = 0.03;           % relative change in foot arch length at mtp 20� dorsiflexion

%% Midtarsal joint
% This will always have the windlass mechanism.
% To simulate a case without windlass, set PF stiffness to 'none' and make
% the torsion spring representing the other ligaments sufficiently stiff,
% also set mtp to spring.

S.mtj = 1;              % 1: use a model with tmt joint (will override tmt)
% plantar fascia
S.PF_stiffness = 'Natali2010'; % stiffness model for the gait simulation
        % options:
        % 'none''linear''Gefen2002''Cheng2008''Natali2010''Song2011'
S.sf_PF = 10;                % multiply PF force with constant factor
S.PF_slack_length = 0.15; % (m) slack length

S.R_mtth = 9.5e-3;
S.WLpoly = 1;

% other ligaments (long, short planter ligament, etc)
S.MT_li_nonl = 0;       % 1: nonlinear torque-angle characteristic
% S.mtj_stiffness = 'Gefen2002';
% S.mtj_stiffness = 'Ker1987';
% S.mtj_stiffness = 'signed_lin';
S.mtj_stiffness = 'fitted6';

S.kMT_li = 300;          % angular stiffness in case of linear
S.kMT_li2 = 10;          % angular stiffness in case of linear
% S.dMT = 5;               % (Nms/rad) damping

S.stiffen_arch = 0;      % (Nm/rad) extra stiffness added to arch (mtj)

% PF reaction torque on mtp joint
S.WL_T_mtp = 1;         % 0: spring mtp, 1: PF reaction on mtp
S.Mu_mtp = 1;           % 0: torque actuator, 1: muscles connected to mtp
    
S.kMTP = 1;
% S.dMTP = 0;

% List of stiffness models to use for the STATIC footmodel:
PF_stiffness = {S.PF_stiffness};
% PF_stiffness = {'linear','Gefen2002','Cheng2008','Barrett2018','Natali2010','none'};
% PF_stiffness = {'none'};


%% Exoskeleton
% exo
S.ExoBool       = 0;    % 1: is wearing exo
S.ExoScale      = 0;    % scale factor of exoskeleton assistance profile 
                        % 0: no assistance (passive) 1: nominal assistance (active)
                        
S.DataSet = 'PoggenSee2020_AFO';            % dataset with exoskeleton torque profile

% Where the exo torque is transferred to the ankle
S.ExoImplementation = 'TorqueTibiaCalcn';
% S.ExoImplementation = 'TorqueTibiaCalcnMetatarsi';
% S.ExoImplementation = 'TorqueTibiaMetatarsi';

% Ideal assistance
ia = 0;
% S.T_max_ankle_exo = 30;
% S.T_min_ankle_exo = 0;
% S.P_max_ankle_exo = 50;



%% Initial guess
% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 3;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)

if S.IGmodeID == 4
    S.savename_ig   = 'NoExo';
elseif S.IGmodeID == 3
    S.ResultsF_ig   = 'Final';
    if strcmp(S.subject,'s1_Poggensee')
        S.savename_ig   = 'Pog_s1_bCst_ig24';
    else
        S.savename_ig   = 'Fal_s1_bCst_ig21';
%         S.savename_ig   = 'Fal_s1_bCst_ig1_v27';
    end
end

if S.IGmodeID == 1 || S.IGsel
    if strcmp(S.subject,'s1_Poggensee')
        S.IG_PelvisY = 0.896 + 0.0131;
    else
        S.IG_PelvisY = 0.9385 + 0.0131;
    end
end

%% Automated settings
% Change some more settings, based on what was selected above. 
% No need to change anything below this line to run the code.

% Ideal assistance
if ia
    S.ExoController = 'Ideal Assistance';
end

% select folder with polynomials
S.PolyFolder = S.subject;

% external function
if S.tmt == 0 && S.mtj == 0
    if strcmp(S.subject,'s1_Poggensee')
        if S.ExoBool == 0
            S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtp.dll';
            S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtp_pp.dll';
        else
            S.ExternalFunc  = 'SimExo_3D_talus_out.dll';
        end
    elseif strcmp(S.subject,'subject1')
        if S.ExoBool == 0
%             S.ExternalFunc  = 'ID_Subject1.dll';
%             S.ExternalFunc2 = 'PredSim_3D_Fal_s1_pp_v2.dll';
            S.ExternalFunc  = 'PredSim_3D_Fal_s1_v7.dll';
            S.ExternalFunc2 = 'PredSim_3D_Fal_s1_pp_v6.dll';
        end
    end
    
elseif S.mtj == 1
    if S.ExoBool == 0
        if strcmp(S.subject,'s1_Poggensee')
            S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtj_v3.dll';
            S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtj_pp_v3.dll';
            
        elseif strcmp(S.subject,'subject1')
             if S.contactStiff == 10
                S.ExternalFunc  = 'PredSim_3D_Fal_s1_mtj_spx10_v1.dll';
                S.ExternalFunc2 = 'PredSim_3D_Fal_s1_mtj_spx10_pp_v3.dll';
             elseif S.contactStiff == 2
                S.ExternalFunc  = 'PredSim_3D_Fal_s1_mtj_spx2_v1.dll';
                S.ExternalFunc2 = 'PredSim_3D_Fal_s1_mtj_spx2_pp_v2.dll';
             else
                S.ExternalFunc  = 'PredSim_3D_Fal_s1_mtj_v1.dll';
                S.ExternalFunc2 = 'PredSim_3D_Fal_s1_mtj_pp_v6.dll';
             end
        end
    end
    
elseif S.tmt == 1
    if S.ExoBool == 0
        if strcmp(S.subject,'s1_Poggensee')
            S.ExternalFunc  = 'PredSim_3D_Pog_s1_tmt_v3.dll';
            S.ExternalFunc2 = 'PredSim_3D_Pog_s1_tmt_pp_v3.dll';
            
        elseif strcmp(S.subject,'subject1')
            S.ExternalFunc  = 'PredSim_3D_Fal_s1_tmt_v1.dll';
            S.ExternalFunc2 = 'PredSim_3D_Fal_s1_tmt_pp_v1.dll';
        end
    else
        if strcmp(S.ExoImplementation,'TorqueTibiaCalcn')            
            S.ExternalFunc  = 'SimExo_3D_Pog_s1_tmt_TTC_v3.dll';
            S.ExternalFunc2  = 'SimExo_3D_Pog_s1_tmt_TTC_pp_v3.dll';
            
        elseif strcmp(S.ExoImplementation,'TorqueTibiaCalcnMetatarsi')
            S.ExternalFunc  = 'SimExo_3D_Pog_s1_tmt_TTCM_v1.dll';
            S.ExternalFunc2  = 'SimExo_3D_Pog_s1_tmt_TTCM_pp_v1.dll';
            
        elseif strcmp(S.ExoImplementation,'TorqueTibiaMetatarsi')
            S.ExternalFunc  = 'SimExo_3D_Pog_s1_tmt_TTM_v1.dll';
            S.ExternalFunc2  = 'SimExo_3D_Pog_s1_tmt_TTM_pp_v1.dll';
        end
    end
    
end

% build standardised names
[savename, casfuncfol] = getSavename(S);
S.CasadiFunc_Folders = [casfuncfol suffixCasName];
S.savename = [savename suffixName];


if batchQueue
    if (exist([pathRepo '/Results/batchQ.mat'],'file')==2) 
        load([pathRepo '/Results/batchQ.mat'],'batchQ');
    else
        batchQ.(S.savename) = struct('S',[]);
    end
    batchQ.(S.savename).S = S;
    
else
    % make folder to store results if it doesn't exist
    pathResults = fullfile([pathRepo '/Results'],S.ResultsFolder);
    if ~isfolder(pathResults)
        mkdir(pathResults);
    end

    % Create the casadifunctions if they do not exist yet
    if ~isfolder([pathRepo '\CasADiFunctions\' S.CasadiFunc_Folders]) && slv
        disp('Creating casadifunctions...');
        CreateCasADiFunctions_all_tmt(pathRepo,S);
        disp('...casadifunctions created');
    end

end

%% Run

if S.tmt == 0 && S.mtj == 0
    if batchQueue
        batchQ.(S.savename).PredSim = 'f_PredSim_Gait92';
        batchQ.(S.savename).LoadSim = 'f_LoadSim_Gait92';
    end
    if slv        % run the optimization
        f_PredSim_Gait92(S);
    end
    if pp           % post-proces simulation results
        f_LoadSim_Gait92(S.ResultsFolder,S.savename);
    end
else
    if ia
        if batchQueue
            batchQ.(S.savename).PredSim = 'f_PredSim_Gait92_tmt_ia';
            batchQ.(S.savename).LoadSim = 'f_LoadSim_Gait92_tmt_ia';
        end
        if slv
            f_PredSim_Gait92_tmt_ia(S);
        end
        if pp
            f_LoadSim_Gait92_tmt_ia(S.ResultsFolder,S.savename);
        end
    else
        if batchQueue
            if S.Mu_mtp
                batchQ.(S.savename).PredSim = 'f_PredSim_Gait92_tmt_v2';
                batchQ.(S.savename).LoadSim = 'f_LoadSim_Gait92_tmt_v2';
            else
                batchQ.(S.savename).PredSim = 'f_PredSim_Gait92_tmt';
                batchQ.(S.savename).LoadSim = 'f_LoadSim_Gait92_tmt';
            end
        end
        if slv
            if S.Mu_mtp
                f_PredSim_Gait92_tmt_v2(S);
            else
                f_PredSim_Gait92_tmt(S);
            end
        end
        if pp
            if S.Mu_mtp
                f_LoadSim_Gait92_tmt_v2(S.ResultsFolder,S.savename);
            else
                f_LoadSim_Gait92_tmt(S.ResultsFolder,S.savename);
            end
        end
    end
end

if batchQueue
    save([pathRepo '/Results/batchQ.mat'],'batchQ');
end

%% Plot
if plot
    h = figure();
    set(h,'Position',[82 151 1497 827]);
    PlotResults_3DSim_tmt(fullfile(pathResults,[S.savename '_pp.mat']),[1 0 0],savename,h);
end


%% Static model of foot
if foot_standing
    for i=1:numel(PF_stiffness)
        S.PF_stiffness = PF_stiffness{i};
%         f_staticFootCompression_v2(S);
        f_staticFootCompression_v5(S);
    end
end
if foot_hanging
    for i=1:numel(PF_stiffness)
        S.PF_stiffness = PF_stiffness{i};
%         f_staticFootHanging(S);
        f_staticFootHanging_v2(S);
    end
end



