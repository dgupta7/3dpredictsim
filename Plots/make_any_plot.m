clear all
clc

%% Paths
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/VariousFunctions']);

%% Settings
% Folder will be filtered to only plot results that satisfy all chosen
% settings. Put an entry in comment to not use it to filter.

% folder to filter from
ResultsFolder = 'Batchsim_tmt_linear_v2'; %'Batchsim_tmt_linear_v2' 'Batchsim_tmt_linear' 'all'

% experimental data to plot as reference
reference_data = 'norm'; % 'none' 'norm' 'pas' act'

% tarsometatarsal joint
S.tmt = 1;              % 1: use a model with tmt joint
S.tmt_locked = 0;       % 1: lock the tmt joint (to compare with model w/o)
% S.kTMT = 800;           % [250 500 800 1000 2000] (Nm/rad) stiffness of tmt joint 
S.dTMT = 0.2;             % [0 0.2 0.5] (Nms/rad) damping of tmt joint


% assumption to simplify Hill-type muscle model
S.MuscModelAsmp = 0;    % 0: musc height = cst, 1: pennation angle = cst

% exo
S.ExoBool       = 0;    % 1: is wearing exo
S.ExoScale      = 0;    % scale factor of exoskeleton assistance profile 
                        % 0: no assistance (passive) 1: nominal assistance (active)
% S.IGsel         = 2;    % initial guess identifier (1: quasi random, 2: data-based)
% S.IGmodeID      = 4;    % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)



pathResult = fullfile([pathRepo '/Results/' ResultsFolder]);
pathData = [pathRepo,'/ExperimentalData','/ExperimentalData.mat'];



% build criteria to filter resultsfolder
ct=1;
if isfield(S,'subject') && ~isempty(S.subject) && strcmp(S.subject,'s1_Poggensee')
    criteria{ct} = 'Pog_s1';
    ct=ct+1;
end
if isfield(S,'tmt') && ~isempty(S.tmt)
    if S.tmt && S.tmt_locked
        criteria{ct} ='tmtL';
        ct=ct+1;
    end
    if S.tmt && ~S.tmt_locked
        criteria{ct} ='tmt';
        ct=ct+1;
        criteria{ct} ='not_tmtL';
        ct=ct+1;
    end
    if ~S.tmt
        criteria{ct} ='not_tmt';
        ct=ct+1;
    end
end
if isfield(S,'MuscModelAsmp') && ~isempty(S.MuscModelAsmp) 
    if S.MuscModelAsmp==0
    criteria{ct} = 'bCst';
    ct=ct+1;
    elseif S.MuscModelAsmp==1
        criteria{ct} = 'aCst';
    ct=ct+1;
    end
end
if isfield(S,'dTMT') && ~isempty(S.dTMT) && ~S.tmt_locked
    criteria{ct} = ['d0' num2str(S.dTMT*10)];
    ct=ct+1;
end
if isfield(S,'kTMT') && ~isempty(S.kTMT) && ~S.tmt_locked
    criteria{ct} = ['k' num2str(S.kTMT)];
    ct=ct+1;
end
if isfield(S,'IGsel') && ~isempty(S.IGsel)
    if S.IGsel == 1
        criteria{ct} = 'ig1';
        ct=ct+1;
    else
        criteria{ct} = ['ig2' num2str(S.IGmodeID )];
        ct=ct+1;
    end
end
if S.ExoBool == 1
    if S.ExoScale == 0
        criteria{ct} = 'pas';
        ct=ct+1;
        reference_data = 'pas';
    else
        criteria{ct} = 'act';
        ct=ct+1;
        reference_data = 'act';
    end
else
    criteria{ct} = 'not_pas';
    ct=ct+1;
    criteria{ct} = 'not_act';
    ct=ct+1;
end


%%

[filteredResults] = filterResultfolderByParameters(pathResult,criteria);
n = length(filteredResults);
filteredResults{n+1} = 'D:\school\WTK\thesis\model\3dpredictsim\Results\all\Pog_s1_bCst_pp.mat';

Plot3D(filteredResults,reference_data)

if length(filteredResults)>1
    ValidationPlots(pathData,filteredResults{1},filteredResults{2:end})
else
    ValidationPlots(pathData,filteredResults{1})
end

%%


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