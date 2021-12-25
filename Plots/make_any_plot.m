%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script allows to plot any of the results from the Thesis_Lars
% branch. See \3dpredictsim\Plots\make_preselected_plots.m to generate
% figures preselected groups of results.
% 
% Select one or more folders with results to choose from. Then set
% values to the parameters you want to filter out. Put the paramaters in
% comment if you want to see results for any value. 
% The filter method relies on the name of the resultfile, so this does
% limit the possible filter criteria. 
% Line 127-133 shows how to manually add more filter criteria. It will only
% return filenames which contain that string. Use 'not_string' to get the
% filenames that do not contain 'string'.
% If no file in a folder contains a given string, that criterium will be
% dropped, for that folder only.
%
% Author: Lars D'Hondt (Dec 2021)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
clc

%% Paths
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/VariousFunctions']);
addpath([pathRepo '/PassiveMoments']);
addpath([pathRepo '/FootModel']);



%% Folder(s) with results
ResultsFolder = {'debug'};

%% General information

% S.suffixName = 'v4';

S.subject = 'Fal_s1';
S.Foot.Model = 'mtj';
S.Foot.Scaling = 'custom';

% S.Foot.kMTP = 17;

% S.Foot.MT_li_nonl = 1;       % 1: nonlinear torque-angle characteristic
% S.Foot.mtj_stiffness = 'MG_exp_table';

% S.Foot.PIM = 1;





% initial guess
% S.IGsel         = 2;    % initial guess identifier (1: quasi random, 2: data-based)
% S.IGmodeID      = 3;    % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% results = {
% %     '\debug\Fal_s1_mtj_sc_cspx10_oy_MTPm_k1_d01_MTJm_nl_MG_exp_table_d01_PF_Gefen2002_ls146_PIM_w1e+03_w1e+04_ig21'
% %     '\debug\Fal_s1_mtj_sc_cspx10_oy_MTPm_k1_d01_MTJm_nl_MG_exp_table_d01_PF_Gefen2002_ls146_PIM_w1e+03_w1e+05_ig21'
%     '\debug\Fal_s1_mtj_sc_cspx10_oy_MTPm_k1_d01_MTJm_nl_MG_exp_table_d01_PF_Gefen2002_ls146_PIM_w1e+04_w1e+05_ig21'
%     '\debug\Fal_s1_mtj_sc_cspx10_oy_MTPm_k1_d01_MTJm_nl_MG_exp_table_d01_PF_Gefen2002_ls146_PIM_w1e+05_w1e+05_ig21'
% %     '\debug\Fal_s1_mtj_sc_cspx10_oy_MTPm_k1_d01_MTJm_nl_MG_exp_table_d01_PF_Gefen2002_ls146_PIM_w2e+04_w1e+05_ig21'
%     '\debug\Fal_s1_mtj_sc_cspx10_oy_MTPm_k1_d01_MTJm_nl_MG_exp_table_d01_PF_Gefen2002_ls146_PIM2_w1e+05_w1e+05_ig21'
%     };


if exist('results','var') && ~isempty(results)
    filteredResults = {length(results)};
    for i=1:length(results)
        filteredResults{i} = fullfile(pathRepo, 'Results', [results{i} '_pp.mat']);
    end
else

    % get file names
    pathResult = {numel(ResultsFolder)};
    for i=1:numel(ResultsFolder)
        pathResult{i} = fullfile([pathRepo '/Results/' ResultsFolder{i}]);
    end
    
    
    % get filter criteria
    [~,~,criteria] = getSavename(S);
    
    criteria{end+1} = 'cspx10_oy';
    criteria{end+1} = 'not_PIM';
    
    
    % filter filenames
    [filteredResults] = filterResultfolderByParameters(pathResult,criteria);
    
    for i=1:numel(filteredResults)
        rfpathi = filteredResults{i};
        idx_i = strfind(rfpathi,'\');
        disp(rfpathi(idx_i(end-1):end));
    end

end
ref = {};

% ref{end+1} = fullfile([pathRepo '/Results/debug\Fal_s1_mtp_sd_MTPp_k17_d05_ig21_pp.mat']);
% ref{end+1} = fullfile([pathRepo '/Results/debug\Fal_s1_mtp_sc_MTPp_k17_d05_ig21_pp.mat']);
% ref{end+1} = fullfile([pathRepo '/Results/debug\Fal_s1_mtp_sd_MTPm_k1_d05_ig21_pp.mat']);


% filteredResultsWithRef = filteredResults';
% filteredResultsWithRef = [ref, filteredResults]';
filteredResultsWithRef = [filteredResults, ref]';


ResultsFile = filteredResultsWithRef;

LegNames = {''};

mtj = 1;
figNamePrefix = 'none';

%%% select figures to make
makeplot.kinematics                     = 1; % selected joint angles
makeplot.kinetics                       = 0; % selected joint torques
makeplot.ankle_musc                     = 1; % ankle muscles
makeplot.GRF                            = 0; % ground interaction
makeplot.compareLiterature              = 0; % mtj and mtp Caravaggi 2018
makeplot.compareTakahashi17             = 0; % "distal to segment" power analysis
makeplot.compareTakahashi17_separate    = 0; % "distal to segment" power analysis
makeplot.compareTakahashi17_mtj_only    = 0; % plot mtj power over experimental result
makeplot.compareTakahashi17_W_bar       = 0; % "distal to segment" work analysis
makeplot.allQsTs                        = 0; % all joint angles and torques
makeplot.windlass                       = 1; % plantar fascia and foot arch info
makeplot.power                          = 0; % datailed power decomposition
makeplot.work                           = 0; % same as power, but work over GC
makeplot.work_bar                       = 0; % positive, negative and net work bar plot
makeplot.work_bar_small                 = 0; % positive, negative and net work bar plot
makeplot.power_main                     = 0; % main power components of foot
makeplot.spatiotemp                     = 0; % stridelength etc.
makeplot.ankle_correlation              = 0; % correlation of ankle 
makeplot.E_muscle_bar                   = 0; % muscle metabolic energy totals
makeplot.E_muscle_bar_small             = 0; % metabolic energy and work by selected muscle groups
makeplot.toes                           = 0; % toe flexor and extensor muscle info
makeplot.Edot_all                       = 0; % summed metabolic energy rate
makeplot.Energy_cost                    = 0; % decompose metabolic cost components
makeplot.muscle_act                     = 0; % muscle activity
makeplot.muscle_joint_moment            = 1; % moments of muscles around ankle-foot joints
makeplot.muscle_joint_power             = 0; % powers of muscles around ankle-foot joints
makeplot.Objective_cost                 = 0; % cost function decomposition

PlotResults_3DSim_Report(ResultsFile,LegNames,'Fal_s1',mtj,makeplot,figNamePrefix);


