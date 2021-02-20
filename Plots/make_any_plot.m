clear all
clc

%% Paths
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/VariousFunctions']);
addpath([pathRepo '/PassiveMoments']);

%% Settings
% Folder will be filtered to only plot results that satisfy all chosen
% settings. Put an entry in comment to not use it to filter.

plot_default = 1;
plot_validation = 0;

% folder to filter from
ResultsFolder = 'debug'; % 'tmt_lin' 'debug' 'debug_batch' 'running'

% experimental data to plot as reference
reference_data = 'norm'; % 'none' 'norm' 'pas' 'act' 'Fal_s1'


% tarsometatarsal joint
S.tmt = 1;              % 1: use a model with tmt joint
S.tmt_locked = 0;       % 1: lock the tmt joint (to compare with model w/o)
S.kTMT = 1000;           % [250 500 800 1000 2000] (Nm/rad) stiffness of tmt joint 
S.dTMT = 0;             % [0 0.2 0.5] (Nms/rad) damping of tmt joint

S.Windlass = 0;
% S.cWL = 0.02;           % relative change in foot arch length at mtp 20° dorsiflexion

% assumption to simplify Hill-type muscle model
S.MuscModelAsmp = 0;    % 0: musc height = cst, 1: pennation angle = cst

% exo
S.ExoBool       = 1;    % 1: is wearing exo
S.ExoScale      = 1;    % scale factor of exoskeleton assistance profile 
                        % 0: no assistance (passive) 1: nominal assistance (active)

% initial guess
% S.IGsel         = 2;    % initial guess identifier (1: quasi random, 2: data-based)
% S.IGmodeID      = 4;    % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)



pathResult = fullfile([pathRepo '/Results/' ResultsFolder]);
pathData = [pathRepo,'/ExperimentalData','/ExperimentalData.mat'];


if isfield(S,'ExoBool') && ~isempty(S.ExoBool) && S.ExoBool == 1
    if S.ExoScale == 0 && ~strcmp(reference_data,'none')
        reference_data = 'pas';
    elseif ~strcmp(reference_data,'none')
        reference_data = 'act';
    end
end


[~,~,criteria] = getSavename(S);

criteria{end+1} = 'ia';

%%

[filteredResults] = filterResultfolderByParameters(pathResult,criteria);

% use this as reference
n = length(filteredResults);
if isfield(S,'ExoBool') && ~isempty(S.ExoBool) && S.ExoBool == 1
    if S.ExoScale == 1
        ref{1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_act_pp.mat';
    else
        ref{1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_pas_pp.mat';
    end
else
    ref{1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_pp.mat';
end

% ref{2} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\tmt_lin\Pog_s1_tmtL_bCst_ig24_v3_pp.mat';
% ref{2} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\debug_batch\Pog_s1_tmt_bCst_d00_k1000_ig24_pp.mat';

% ref = {'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_pp.mat',...
%        'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_pas_pp.mat',...
%        'D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst_act_pp.mat'};

% filteredResults{n+1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\debug_tmt\Pog_s1_tmt_bCst_d02_k800_kc1_t5_ig24_v3_pp.mat';
% filteredResults{n+2} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\debug_tmt\Pog_s1_tmt_bCst_d02_k800_ig24_v3_pp.mat';

filteredResultsWithRef = {ref{:}, filteredResults{:}};

%%

if plot_default
%     Plot3D(filteredResultsWithRef,reference_data)
    Plot3D(filteredResults,reference_data)
%     Plot3D(ref,reference_data)
end

if plot_validation
    if strcmp(reference_data,'Fal_s1')
        pathData = 'D:\school\WTK\thesis\model\3dpredictsim\Data\Fal_s1.mat';
    else
        pathData = 'D:\school\WTK\thesis\model\3dpredictsim\Data\Pog_s1.mat';
    end
    if length(filteredResultsWithRef)>1
        ValidationPlots(pathData,filteredResultsWithRef{1},filteredResultsWithRef{2:end})
    else
        ValidationPlots(pathData,filteredResultsWithRef{1})
    end
end


%% old/temp stuff


% Plot3D_pwd_separate(pathResult); % plot default figures for entire resultsfolder

%%
% if plot_comp
%     ResultsFile1 = fullfile(pathResults,'Pog_s1_tmt_bCst_d02_k800_ig24_pp.mat');
%     ResultsFile2 = fullfile(pathResults,'Pog_s1_tmt_bCst_d02_k1000_ig24_pp.mat');
% 
%     % plot default figures
%     hh1 = figure();
%     set(hh1,'Position',[82         151        1497         827]);
%     PlotResults_3DSim_tmt(ResultsFile1,[1 0 0],'k800',hh1);
%     PlotResults_3DSim_tmt(ResultsFile2,[0 0 1],'k1000',hh1);
%     
%     % Compare kinematics and kinetics any number of results to the measurementdata
%     ValidationPlots(pathData,ResultsFile1,ResultsFile2);
% 
%     % Compare muscles for 2 results
% %     PlotResultsComparison_3DSim(ResultsFile1,ResultsFile2);
% 
% end
% 
% if timing
%     ResultsFile = fullfile(pathResults,'Pog_s1_tmt_bCst_d02_k800_ig24_pp.mat');
%     Plot_tmt_mtp(ResultsFile)
% end