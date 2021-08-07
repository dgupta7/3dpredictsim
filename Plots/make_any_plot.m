%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script allows to plot any of the results from the Thesis_Lars
% branch. Select one or more folders with results to choose from. Then set
% values to the parameters you want to filter out. Put the paramaters in
% comment if you want to see results for any value. 
% The filter method relies on the name of the resultfile, so this does
% limit the possible filter criteria. 
% Line 127-133 shows how to manually add more filter criteria. It will only
% return filenames which contain that string. Use 'not_string' to get the
% filenames that do not contain 'string'.
% If no file in a folder contains a given string, that criterium will be
% dropped, for that file only.
%
% Author: Lars D'Hondt (May 2021)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%% Paths
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/VariousFunctions']);
addpath([pathRepo '/PassiveMoments']);
addpath([pathRepo '/FootModel']);

%% Settings
plot_default = 0;
plot_validation = 0;
plot_report = 0;

plot_foot_standing = 1;
plot_foot_hanging = 0;


%% Folder(s) with results
% ResultsFolder = {'batch_windlass','batch_tmt_lin'}; % 'tmt_lin' 'debug' 'debug_batch' 'running'
% ResultsFolder = {'batch_windlass'};
% ResultsFolder = {'running'};
% ResultsFolder = {'MuscleModel'};
% ResultsFolder = {'batch_tmt_lin'};
% ResultsFolder = {'test_WL_v2'};
ResultsFolder = {'MidTarsalJoint'};
% ResultsFolder = {'Final'};

%% General information
% experimental data to plot as reference
reference_data = 'norm'; % 'none' 'norm' 'pas' 'act' 'Fal_s1'

% assumption to simplify Hill-type muscle model
S.MuscModelAsmp = 0;    % 0: musc height = cst, 1: pennation angle = cst

% Test subject
S.subject            = 'subject1';
% S.subject            = 's1_Poggensee';

% initial guess
% S.IGsel         = 2;    % initial guess identifier (1: quasi random, 2: data-based)
% S.IGmodeID      = 3;    % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)

%% Footmodel only
% subtalar orientation (footmodel only)
% S.subt_orientation = 'default'; %'Reule2010''Parr2012'

% max compressive force
Fmax = 960;

%% tarsometatarsal joint
% S.tmt = 0;              % 1: use a model with tmt joint
% S.tmt_locked = 0;       % 1: lock the tmt joint (to compare with model w/o)
% % S.kTMT = 1000;          % [250 500 800 1000 2000] (Nm/rad) stiffness of tmt joint 
% % S.dTMT = 0.5;           % [0 0.2 0.5] (Nms/rad) damping of tmt joint
% % Windlass mechanism
% S.Windlass = 1;         % 1: has windlass mechnism
% S.cWL = 0.03;           % relative change in foot arch length at mtp 20� dorsiflexion


%% Midtarsal joint
% This will always have the windlass mechanism.
S.mtj = 1;              % 1: use a model with tmt joint (will override tmt)
% plantar fascia
% S.PF_stiffness = 'Song2011'; % stiffness model for the gait simulation
% S.PF_stiffness = 'Natali2010';
% S.PF_stiffness = 'None';
        % options: 'none''linear''Gefen2002''Cheng2008''Barrett2018''Natali2010''Song2011'
S.PF_slack_length = 0.150; % slack length (m)
% other ligaments (long, short planter ligament, etc)
% S.MT_li_nonl = 1;       % 1: nonlinear torque-angle characteristic
% S.kMT_li = 300;         % angular stiffness in case of linear
% S.mtj_stiffness = 'Gefen2001';
% S.mtj_stiffness = 'Ker1987';
% S.mtj_stiffness = 'Song2011';
% S.mtj_stiffness = 'signed_lin';
% S.dMT = 0;

% PF reaction torque on mtp joint
S.WL_T_mtp = 1;         % 0: spring mtp, 1: PF reaction on mtp
S.Mu_mtp = 0;           % 0: torque actuator, 1: muscles connected to mtp
% S.kMTP = 5;
S.dMTP = 0;

%% Exoskeleton
S.ExoBool       = 0;    % 1: is wearing exo
S.ExoScale      = 0;    % scale factor of exoskeleton assistance profile 
                        % 0: no assistance (passive) 1: nominal assistance (active)
% S.ExoImplementation = 'TorqueTibiaCalcn';
% S.ExoController = 'Ideal Assistance';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if plot_validation || plot_default || plot_report
    % get file names
    for i=1:numel(ResultsFolder)
        pathResult{i} = fullfile([pathRepo '/Results/' ResultsFolder{i}]);
    end
    
    pathData = [pathRepo,'/ExperimentalData','/ExperimentalData.mat'];

    % set reference data
    if isfield(S,'ExoBool') && ~isempty(S.ExoBool) && S.ExoBool == 1
        if S.ExoScale == 0 && ~strcmp(reference_data,'none')
            reference_data = 'pas';
        elseif ~strcmp(reference_data,'none')
            reference_data = 'act';
        end
    end

    % get filter criteria
    [~,~,criteria] = getSavename(S);

    % manually add more filter criteria
%     criteria{end+1} = 'v27';
%     criteria{end+1} = 'not_Mu1_ig24';
%     criteria{end+1} = 'not_k10';
%     criteria{end+1} = 'spx';
%     criteria{end+1} = 'PFx';
%     criteria{end+1} = 'not_PFx2';
%     criteria{end+1} = 'not_d0';

    % filter filenames
    [filteredResults] = filterResultfolderByParameters(pathResult,criteria);

    for i=1:numel(filteredResults)
        disp(filteredResults{i});
    end
    

    % specify reference results
    n = length(filteredResults);
    if isfield(S,'ExoBool') && ~isempty(S.ExoBool) && S.ExoBool == 1
        if S.ExoScale == 1
            ref{1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_act_pp.mat';
        else
            ref{1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_pas_pp.mat';
        end
    else
        if isfield(S,'subject') && strcmp(S.subject,'subject1')
%             ref{1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Fal_s1_bCst_ig24_v2_pp.mat';
%             ref{1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\test_pp\Fal_s1_bCst_ig24_v2_pp.mat';
            ref{1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\Final\Fal_s1_bCst_ig1_pp.mat';
        elseif isfield(S,'subject') && strcmp(S.subject,'s1_Poggensee')
            ref{1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_pp.mat';
        else
            ref{1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Fal_s1_bCst_ig24_v2_pp.mat';
            ref{2} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_pp.mat';
        end
    end



    filteredResultsWithRef = {ref{:}, filteredResults{:}};
%     filteredResultsWithRef = {filteredResults{:}, ref{:}};

    % set value according to which figure(s) to make
    pl = 0;
    if plot_validation
        if plot_default
            pl = 2;
        else
            pl = 1;
        end
    end
    
    if plot_validation || plot_default
    % call plot function
        [h_default,h_validation] = Plot3D(filteredResultsWithRef,reference_data,pl);
%         Plot3D(filteredResults,reference_data,pl)
%         Plot3D(ref,reference_data,pl)
    end
    
end

%%
if plot_report

%     ResultsFile = {'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Fal_s1_bCst_ig24_v2_pp.mat'};
%         'D:\school\WTK\thesis\model\3dpredictsim\Results\batch_windlass\Fal_s1_tmt_bCst_d05_k1000_WL30_ig24_pp.mat'};
    
    ResultsFile = filteredResultsWithRef;
%     LegNames = {'original','k = 30, PF x10','k = 50, PF x5','k = 50'};
%     LegNames = {'no mtj','50','100','150','200','250','300','400','500','800'};
%     LegNames = {'original','Gefen2002','Natali2010','Song2011','linear'};
%     ResultsFile = filteredResults;
%     LegNames = {'mtj and wl'};
    LegNames = {'mtp', 'mtp, mtj and windlass'};
    
%     ResultsFile = ref;
%     LegNames = {'Falisse 2019'};

    RefData = 'Fal_s1';
    mtj = 1;
    
    makeplot.kinematics = 0;
    makeplot.kinetics = 0;
    makeplot.soleus = 0;
    makeplot.GRF = 0;
    makeplot.compareLiterature = 0;
    makeplot.COP = 0;
    makeplot.allQsTs = 0;
    makeplot.k_mtj_lin = 0;
    makeplot.windlass = 0;
    makeplot.power = 1;
    makeplot.power_T = 0;

    figNamePrefix = 'none';
%     figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\matlab\conclusion';
%     figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\matlab_final\SOTA';

    
    PlotResults_3DSim_Report(ResultsFile,LegNames,RefData,mtj,makeplot,figNamePrefix);
end



%% Static foot model


if plot_foot_standing || plot_foot_hanging
    % get filter criteria
    crit = {};
    if isfield(S,'subject')
        if strcmp(S.subject,'s1_Poggensee')
            crit{end+1} = 'Pog_s1';
        elseif strcmp(S.subject,'subject1')
            crit{end+1} = 'Fal_s1';
        end
    end
    if isfield(S,'subt_orientation')
        if strcmp(S.subt_orientation,'Reule2010')
            crit{end+1} = 'subt2';
        elseif strcmp(S.subt_orientation,'Parr2012')
            crit{end+1} = 'subt3';
        else
            crit{end+1} = 'subt1';
        end
    end
    if plot_foot_standing
        crit{end+1} = 'not_hanging';
    else
        crit{end+1} = 'hanging';
    end
    if isfield(S,'PF_stiffness')
        crit{end+1} = S.PF_stiffness;
    end
    if exist('Fmax','var')
        crit{end+1} = ['_' num2str(Fmax)];
    end
    
%     crit{end+1} = '_sb';
%     crit{end+1} = '_ls150';
%     crit{end+1} = '_Q0_30';
    
    OutFolder{1} = fullfile(pathRepo,'Results','FootModel');
    
    % filter resultfiles
    [resultFiles] = filterResultfolderByParameters(OutFolder,crit);
    
    for i=1:numel(resultFiles)
        disp(resultFiles{i});
    end
    
%     resultFiles = {'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Gefen2002_Gefen2002_Q-30_30_F0_3000_WLv3_ls150.mat',
%         'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_none_Gefen2002_Q-30_30_F0_3000_WLv3_ls150.mat',
%         'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_none_Ker1987_Q-30_30_F0_3000_WLv3_ls150.mat'};
    
%     resultFiles = {'Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_fitted4_Q-20_30_F0_0_WLv3_ls148_mtp1.mat',
%         'Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_fitted4_Q-20_30_F0_0_WLv3_ls148_mtp2.mat'};
    
%     resultFiles = {'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k50_Q0_30_F0_1000_WLv3_ls150_sb1.mat',
%         'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k300_Q0_30_F0_1000_WLv3_ls150_sb1.mat',
%         'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k800_Q0_30_F0_1000_WLv3_ls150_sb1.mat',
%         'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_linear_k300_Q0_30_F0_1000_WLv3_ls150_sb1.mat',
%         'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Gefen2002_Gefen2002_Q0_30_F0_1000_WLv3_ls150_sb1.mat',
%         'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Song2011_Song2011_Q0_30_F0_1000_WLv3_ls150_sb1.mat'};
    
% resultFiles = {'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k50_Q0_30_F0_960_WLv3_ls150_sb1.mat'
% 'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k100_Q0_30_F0_960_WLv3_ls150_sb1.mat'
% 'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k150_Q0_30_F0_960_WLv3_ls150_sb1.mat'
% 'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k200_Q0_30_F0_960_WLv3_ls150_sb1.mat'
% 'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k250_Q0_30_F0_960_WLv3_ls150_sb1.mat'
% 'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k300_Q0_30_F0_960_WLv3_ls150_sb1.mat'
% 'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k350_Q0_30_F0_960_WLv3_ls150_sb1.mat'
% 'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k380_Q0_30_F0_960_WLv3_ls150_sb1.mat'
% 'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k400_Q0_30_F0_960_WLv3_ls150_sb1.mat'
% 'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k450_Q0_30_F0_960_WLv3_ls150_sb1.mat'
% 'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k500_Q0_30_F0_960_WLv3_ls150_sb1.mat'
% 'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k800_Q0_30_F0_960_WLv3_ls150_sb1.mat'
% 'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k1000_Q0_30_F0_960_WLv3_ls150_sb1.mat'};

resultFiles = {'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k50_Q0_30_F0_960_WLv3_ls150_sb1.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k50_Q0_30_F0_960_WLv3_ls150_sb1.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k100_Q0_30_F0_960_WLv3_ls150_sb1.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k150_Q0_30_F0_960_WLv3_ls150_sb1.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k200_Q0_30_F0_960_WLv3_ls150_sb1.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k250_Q0_30_F0_960_WLv3_ls150_sb1.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k300_Q0_30_F0_960_WLv3_ls150_sb1.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k400_Q0_30_F0_960_WLv3_ls150_sb1.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k500_Q0_30_F0_960_WLv3_ls150_sb1.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k800_Q0_30_F0_960_WLv3_ls150_sb1.mat'};

    % call plot function
    nrf = numel(resultFiles);
    CsV = hsv(nrf);
    for i=1:nrf
        load(resultFiles{i},'R');
        if i==1
            h = PlotResults_FootSim(R,CsV(i,:));
        else
            PlotResults_FootSim(R,CsV(i,:),h);
        end
    end
    
end




