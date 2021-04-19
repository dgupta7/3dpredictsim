clear all
clc

%% Paths
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/VariousFunctions']);
addpath([pathRepo '/PassiveMoments']);
addpath([pathRepo '/FootModel']);

%% Settings
% Folder will be filtered to only plot results that satisfy all chosen
% settings. Put an entry in comment to not use it to filter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_default = 1;
plot_validation = 0;
plot_report = 0;

plot_foot_standing = 0;
plot_foot_hanging = 0;


%% Folder(s) with results
% ResultsFolder = {'batch_windlass','batch_tmt_lin'}; % 'tmt_lin' 'debug' 'debug_batch' 'running'
% ResultsFolder = {'batch_windlass'};
% ResultsFolder = {'running'};
% ResultsFolder = {'MuscleModel'};
% ResultsFolder = {'batch_tmt_lin'};
% ResultsFolder = {'test_WL_v2'};
ResultsFolder = {'MidTarsalJoint'};

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
S.subt_orientation = 'default'; %'Reule2010''Parr2012'

% max compressive force
% Fmax = 4000;

%% tarsometatarsal joint
% S.tmt = 1;              % 1: use a model with tmt joint
% S.tmt_locked = 0;       % 1: lock the tmt joint (to compare with model w/o)
% S.kTMT = 1000;          % [250 500 800 1000 2000] (Nm/rad) stiffness of tmt joint 
% S.dTMT = 0.5;           % [0 0.2 0.5] (Nms/rad) damping of tmt joint
% % Windlass mechanism
% S.Windlass = 1;         % 1: has windlass mechnism
% S.cWL = 0.03;           % relative change in foot arch length at mtp 20° dorsiflexion


%% Midtarsal joint
% This will always have the windlass mechanism.
S.mtj = 1;              % 1: use a model with tmt joint (will override tmt)
% plantar fascia
% S.PF_stiffness = 'Gefen2001'; % stiffness model for the gait simulation
S.PF_stiffness = 'Natali2010';
        % options: 'none''linear''Gefen2001''Cheng2008''Barrett2018''Natali2010'
S.PF_slack_length = 0.148; % slack length (m)
% other ligaments (long, short planter ligament, etc)
S.MT_li_nonl = 1;       % 1: nonlinear torque-angle characteristic
% S.kMT_li = 200;         % angular stiffness in case of linear
% S.mtj_stiffness = 'Gefen2001';
% S.mtj_stiffness = 'Ker1987';
S.mtj_stiffness = 'fitted';

% PF reaction torque on mtp joint
S.WL_T_mtp = 1;         % 0: spring mtp, 1: PF reaction on mtp
S.Mu_mtp = 0;           % 0: torque actuator, 1: muscles connected to mtp


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
    criteria{end+1} = 'not_Wmtp';
%     criteria{end+1} = 'not_kmtp5';

    % filter filenames
    [filteredResults] = filterResultfolderByParameters(pathResult,criteria);

%     for i=1:numel(filteredResults)
%         disp(filteredResults{i});
%     end
    

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
            ref{1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\test_pp\Fal_s1_bCst_ig24_v2_pp.mat';
        elseif isfield(S,'subject') && strcmp(S.subject,'s1_Poggensee')
            ref{1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_pp.mat';
        else
            ref{1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Fal_s1_bCst_ig24_v2_pp.mat';
            ref{2} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_pp.mat';
        end
    end

%     ref{2} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\tmt_lin\Pog_s1_tmtL_bCst_ig24_v3_pp.mat';
%     ref{3} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\batch_tmt_lin\Pog_s1_tmt_bCst_d05_k1000_ig24_pp.mat';
%     ref{2} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\batch_windlass\Pog_s1_tmt_bCst_d05_k1000_WL30_ig24_pp.mat';
%     ref{3} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\batch_windlass\Pog_s1_tmt_bCst_d05_k500_WL30_ig24_pp.mat';
%     ref{4} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\batch_windlass\Pog_s1_tmt_bCst_d05_k500_WL20_ig24_pp.mat';


%     ref = {'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_act_pp.mat',...
%            'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_pas_pp.mat',...
%            'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_pp.mat'};

%     filteredResults{n+1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\debug_tmt\Pog_s1_tmt_bCst_d02_k800_kc1_t5_ig24_v3_pp.mat';
%     filteredResults{n+2} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\debug_tmt\Pog_s1_tmt_bCst_d02_k800_ig24_v3_pp.mat';

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
        Plot3D(filteredResultsWithRef,reference_data,pl)
%         Plot3D(filteredResults,reference_data,pl)
%         Plot3D(ref,reference_data,pl)
    end
    
end

%%
if plot_report

%     ResultsFile = {'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Fal_s1_bCst_ig24_v2_pp.mat'};
%         'D:\school\WTK\thesis\model\3dpredictsim\Results\batch_windlass\Fal_s1_tmt_bCst_d05_k1000_WL30_ig24_pp.mat'};
    
    ResultsFile = filteredResultsWithRef;
%     ResultsFile = filteredResults;
    
%     LegNames = {'Simulated'};
    LegNames = {'without midtarsal joint','midtarsal joint and windlass'};
    
    
    RefData = 'Fal_s1';
    mtj = 1;
    
    makeplot.kinematics = 1;
    makeplot.kinetics = 0;
    makeplot.soleus = 1;
    makeplot.GRF = 0;
    makeplot.compareLiterature = 0;
    makeplot.COP = 1;

    figNamePrefix = 0;%'D:\OneDrive\WTK\thesis\figuren\matlab\SOTA';
    
    PlotResults_3DSim_Report(ResultsFile,LegNames,RefData,mtj,makeplot,figNamePrefix);
end

%% Relative effect COT and stride frequency
% DataFile = 'D:\school\WTK\thesis\model\3dpredictsim\Data\Pog_s1.mat';
% 
% ReferenceFiles = {'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_act_pp.mat',...
%        'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_pas_pp.mat',...
%        'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_pp.mat'};
% 
% 
% ValidationPlots_CotSpatiotemp(filteredResults,ReferenceFiles,DataFile);


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
    
%     crit{end+1} = '_WLv3';
%     crit{end+1} = '_ls';
    
    OutFolder{1} = fullfile(pathRepo,'Results','FootModel');
    
    % filter resultfiles
    [resultFiles] = filterResultfolderByParameters(OutFolder,crit);
    
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




