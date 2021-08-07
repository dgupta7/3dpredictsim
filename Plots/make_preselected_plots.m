clear all
close all
clc

plot_static_foot = 0;
plot_full_separate = 1;
plot_full_tabbed = 0;



%% Paths
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/VariousFunctions']);
addpath([pathRepo '/PassiveMoments']);
addpath([pathRepo '/FootModel']);
ResultsFolder = {'MidTarsalJoint'};

for i=1:numel(ResultsFolder)
    pathResult{i} = fullfile([pathRepo '/Results/' ResultsFolder{i}]);
end

%%
ref{1} = fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig21_pp.mat']);
reference_data = 'norm';

%%

% default, best COT
groupNames{1} = 'd_m_t_j = 0 Nms/rad, k_m_t_p = 5 Nm/rad';
Results_default = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_MTP_T5_ig24_v2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1500_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k2000_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k5000_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k10000_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100000_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1000000_MTP_T5_ig23_pp.mat'])};

% default, all ig
groupNames{end+1} = 'd_m_t_j = 0 Nms/rad, k_m_t_p = 5 Nm/rad (ig21)';
Results_default_ig21 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_MTP_T5_ig21_pp.mat'])};

groupNames{end+1} = 'd_m_t_j = 0 Nms/rad, k_m_t_p = 5 Nm/rad (ig23)';
Results_default_ig23 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k5000_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k10000_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100000_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1000000_MTP_T5_ig23_pp.mat'])};

groupNames{end+1} = 'd_m_t_j = 0 Nms/rad, k_m_t_p = 5 Nm/rad (ig24)';
Results_default_ig24 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_MTP_T5_ig24_v2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1500_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k2000_MTP_T5_ig24_pp.mat'])};


% damping 1, best COT
groupNames{end+1} = 'd_m_t_j = 1 Nms/rad, k_m_t_p = 5 Nm/rad';
Results_d010 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_d010_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_d010_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_d010_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_d010_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_d010_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_d010_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1000_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1500_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k2000_d010_MTP_T5_ig24_pp.mat'])};

% damping 1, all ig
groupNames{end+1} = 'd_m_t_j = 1 Nms/rad, k_m_t_p = 5 Nm/rad (ig21)';
Results_d010_ig21 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_d010_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_d010_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_d010_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_d010_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_d010_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_d010_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_d010_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_d010_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_d010_MTP_T5_ig21_pp.mat'])};

groupNames{end+1} = 'd_m_t_j = 1 Nms/rad, k_m_t_p = 5 Nm/rad (ig23)';
Results_d010_ig23 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_d010_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_d010_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_d010_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_d010_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_d010_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_d010_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_d010_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_d010_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_d010_MTP_T5_ig23_pp.mat'])};

groupNames{end+1} = 'd_m_t_j = 1 Nms/rad, k_m_t_p = 5 Nm/rad (ig24)';
Results_d010_ig24 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1000_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1500_d010_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k2000_d010_MTP_T5_ig24_pp.mat'])};

% damping 2, best COT
groupNames{end+1} = 'd_m_t_j = 2 Nms/rad, k_m_t_p = 5 Nm/rad';
Results_d020 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_d020_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_d020_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_d020_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_d020_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_d020_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_d020_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1000_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1500_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k2000_d020_MTP_T5_ig24_pp.mat'])};

% damping 2, all ig
groupNames{end+1} = 'd_m_t_j = 2 Nms/rad, k_m_t_p = 5 Nm/rad (ig21)';
Results_d020_ig21 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_d020_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_d020_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_d020_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_d020_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_d020_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_d020_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_d020_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_d020_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_d020_MTP_T5_ig21_pp.mat'])};

groupNames{end+1} = 'd_m_t_j = 2 Nms/rad, k_m_t_p = 5 Nm/rad (ig23)';
Results_d020_ig23 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_d020_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_d020_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_d020_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_d020_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_d020_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_d020_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_d020_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_d020_MTP_T5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_d020_MTP_T5_ig23_pp.mat'])};

groupNames{end+1} = 'd_m_t_j = 2 Nms/rad, k_m_t_p = 5 Nm/rad (ig24)';
Results_d020_ig24 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1000_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1500_d020_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k2000_d020_MTP_T5_ig24_pp.mat'])};


groupNames{end+1} = 'd_m_t_j = 5 Nms/rad, k_m_t_p = 5 Nm/rad';
Results_d050 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_d050_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_d050_MTP_T5_ig24_pp.mat'])
%     fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_d050_MTP_T5_ig24_pp.mat']) % not converged
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_d050_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_d050_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_d050_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_d050_MTP_T5_ig24_pp.mat'])};

groupNames{end+1} = 'd_m_t_j = 0 Nms/rad, k_m_t_p = 2 Nm/rad';
Results_T2 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_T2_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_MTP_T2_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_MTP_T2_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_MTP_T2_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_MTP_T2_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_MTP_T2_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T2_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_MTP_T2_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_MTP_T2_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_MTP_T2_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1000_MTP_T2_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1500_MTP_T2_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k2000_MTP_T2_ig24_pp.mat'])};

groupNames{end+1} = 'd_m_t_j = 0 Nms/rad, k_m_t_p = 10 Nm/rad';
Results_T10 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_T10_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_MTP_T10_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_MTP_T10_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_MTP_T10_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_MTP_T10_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_MTP_T10_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T10_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_MTP_T10_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_MTP_T10_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_MTP_T10_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1000_MTP_T10_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1500_MTP_T10_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k2000_MTP_T10_ig24_pp.mat'])};

groupNames{end+1} = 'k_P_F x2';
Results_PFx2 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_T5_ig24_PFx2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_MTP_T5_ig24_PFx2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_MTP_T5_ig24_PFx2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_MTP_T5_ig24_PFx2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_MTP_T5_ig24_PFx2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_MTP_T5_ig24_PFx2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_ig24_PFx2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_MTP_T5_ig24_PFx2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_MTP_T5_ig24_PFx2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_MTP_T5_ig24_PFx2_pp.mat'])
%     fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1000_MTP_T5_ig24_PFx2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1500_MTP_T5_ig24_PFx2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k2000_MTP_T5_ig24_PFx2_pp.mat'])};

groupNames{end+1} = 'vs Falisse 2019';
Results_vs = {
    fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_v8_pp.mat'])
    fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_v10_pp.mat'])
    fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_v12_pp.mat'])
    fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_v14_pp.mat'])
    fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_v16_pp.mat'])
    fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_v18_pp.mat'])
    fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_v27_pp.mat'])};
    
groupNames{end+1} = 'vs Song 2011';
Results_vs_Song = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu0_ig1_v8_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu0_ig1_v10_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu0_ig1_v12_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu0_ig1_v14_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu0_ig1_v16_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu0_ig1_v18_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu0_ig1_v27_pp.mat'])};

groupNames{end+1} = 'PF only';
Results_PFonly = {
%     fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k10_MTP_T5_ig24_PFx10_pp.mat'])
%     fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k10_MTP_T2_ig24_PFx10_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_T5_ig24_PFx10_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_T2_ig24_PFx10_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_MTP_T5_ig24_PFx5_pp.mat'])};

groupNames{end+1} = 'Muscle-driven mtp';
Results_Musc = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu1_ig23_pp.mat'])
%     fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu1_ig24_pp.mat']) % local min
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu1_np_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu1_spx2_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu5_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu5_np_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_Mu2_ig24_PFx10_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_Mu5_ig24_PFx10_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_Mu5_ig24_PFx5_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k20_MTP_Mu5_ig24_PFx5_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_fitted6_MTP_Mu5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu5_ig24_pp.mat'])};
    
groupNames{end+1} = 'k_m_t_j_,_n_e_g = 10 Nm/rad';
Results_singed_lin = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_k100_10_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_k200_10_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_k300_10_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_k400_10_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_k500_10_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_k800_10_MTP_T5_ig24_pp.mat'])};

groupNames{end+1} = 'k_m_t_j_,_n_e_g test';
Results_singed_lin_test = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_T5_ig1_300_50_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_T5_ig24_300_10_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_T5_ig24_300_30_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_T5_ig24_500_30_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_T5_ig24_500_50_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_aCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_T5_ig24_500_50_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_Mu5_ig24_300_10_pp.mat'])};

groupNames{end+1} = 'other PF models';
Results_other_PF = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Gefen2002_ls150_MT_nl_signed_lin_MTP_T5_ig24_300_10_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_signed_lin_MTP_T5_ig24_300_10_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_signed_lin_MTP_T5_ig24_800_50_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_linear_ls150_MT_nl_signed_lin_MTP_T5_ig24_300_10_pp.mat'])};

groupNames{end+1} = 'no PF';
Results_no_PF = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_None_ls150_MT_k250_MTP_T17_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_None_ls150_MT_k500_MTP_T17_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_None_ls150_MT_k800_MTP_T17_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_None_ls150_MT_k1000_MTP_T17_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_None_ls150_MT_k1500_MTP_T17_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_None_ls150_MT_k2000_MTP_T17_ig24_pp.mat'])};
    
groupNames{end+1} = 'd_m_t_p = 0 Nms/rad';
Results_mtp_d0 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_MTP_T5_d00_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_MTP_T5_d00_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_d00_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_MTP_T5_d00_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_MTP_T5_d00_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_MTP_T5_d00_ig24_pp.mat'])};

groupNames{end+1} = 'stiffer contact';
Results_stiff_contact = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_spx10_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu1_np_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu1_spx2_ig24_pp.mat'])};

groupNames{end+1} = 'misc';
Results_misc = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_fitted6_MTP_T2_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_fitted6_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_none_ls148_MT_k1000_MTP_k10_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_none_ls148_MT_k500_MTP_k10_ig24_pp.mat'])
    };

best_COT = fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_ig23_pp.mat']);

total_improvement = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_None_ls150_MT_k1000_MTP_T17_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_T5_ig24_PFx10_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_ig23_pp.mat'])};

PF_stiffening = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_MTP_T5_ig24_v2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_MTP_T5_ig24_PFx2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_MTP_T5_ig24_PFx5_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_T5_ig24_PFx10_pp.mat'])};

%%

filteredResultsWithRef{1}     = {ref{:}, Results_default{:}};           % 1
filteredResultsWithRef{end+1} = {ref{:}, Results_default_ig21{:}};      
filteredResultsWithRef{end+1} = {ref{:}, Results_default_ig23{:}}; 
filteredResultsWithRef{end+1} = {ref{:}, Results_default_ig24{:}}; 


filteredResultsWithRef{end+1} = {ref{:}, Results_d010{:}};              % 5
filteredResultsWithRef{end+1} = {ref{:}, Results_d010_ig21{:}};
filteredResultsWithRef{end+1} = {ref{:}, Results_d010_ig23{:}}; 
filteredResultsWithRef{end+1} = {ref{:}, Results_d010_ig24{:}}; 

filteredResultsWithRef{end+1} = {ref{:}, Results_d020{:}};              % 9
filteredResultsWithRef{end+1} = {ref{:}, Results_d020_ig21{:}};
filteredResultsWithRef{end+1} = {ref{:}, Results_d020_ig23{:}};
filteredResultsWithRef{end+1} = {ref{:}, Results_d020_ig24{:}};

filteredResultsWithRef{end+1} = {ref{:}, Results_d050{:}};              % 13
filteredResultsWithRef{end+1} = {ref{:}, Results_T2{:}};                % 14
filteredResultsWithRef{end+1} = {ref{:}, Results_T10{:}};               % 15
filteredResultsWithRef{end+1} = {ref{:}, Results_PFx2{:}};              % 16
filteredResultsWithRef{end+1} = {ref{:}, Results_vs{:}};                % 17
filteredResultsWithRef{end+1} = {ref{:}, Results_vs_Song{:}};           % 18
filteredResultsWithRef{end+1} = {ref{:}, Results_PFonly{:}};            % 19
filteredResultsWithRef{end+1} = {ref{:}, Results_Musc{:}};              % 20
filteredResultsWithRef{end+1} = {ref{:}, Results_singed_lin{:}};        % 21
filteredResultsWithRef{end+1} = {ref{:}, Results_singed_lin_test{:}};   % 22
filteredResultsWithRef{end+1} = {ref{:}, Results_other_PF{:}};          % 23
filteredResultsWithRef{end+1} = {ref{:}, Results_no_PF{:}};             % 24
filteredResultsWithRef{end+1} = {ref{:}, Results_mtp_d0{:}};            % 25
filteredResultsWithRef{end+1} = {ref{:}, Results_stiff_contact{:}};     % 26
filteredResultsWithRef{end+1} = {ref{:}, Results_misc{:}};              % 27
filteredResultsWithRef{end+1} = {ref{:}, total_improvement{:}};         % 28

idx = [1,5,9,13:16,21]; % for comparison plots in function of k_mtj

%% check for unused results
% Results_all = {};
% for i=1:length(filteredResultsWithRef)
%     Results_all = {Results_all{:},filteredResultsWithRef{i}{2:end}};
% end
% 
% dpath = fullfile([pathRepo '\Results\MidTarsalJoint';
% MatFiles = dir(fullfile(dpath,'*_pp.mat'));
% 
% 
% nFil = length(MatFiles);
% not_yet = ones(nFil,1);
% for i=1:numel(Results_all)
%    for j=1:nFil 
%        if not_yet(j) == 1
%            filename = MatFiles(j).name;
%            if sum(strfind(Results_all{i},filename))>0
%                not_yet(j) = 0;
%            end
%        end
%    end
% end
% 
% for j=1:nFil 
%     if not_yet(j) == 1
%         disp(MatFiles(j).name);
%     end
% end



%% Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Select results to plot

ResultsFile = filteredResultsWithRef{23};
% ResultsFile = {filteredResultsWithRef{20}{[1,7:10]}};

% ResultsFile = {ref{:}};
LegNames = {'Falisse 2019'};

% ResultsFile = {ref{:},best_COT};
% LegNames = {'Falisse 2019','Elastic arch & windlass'};

% ResultsFile = {best_COT};
% LegNames = {'Elastic arch & windlass'};

% ResultsFile = {ref{:},PF_stiffening{:}};

% ResultsFile = filteredResultsWithRef{28};
% LegNames = {'Falisse 2019','Elastic arch','Windlass','Elastic arch & windlass'};


RefData = 'Fal_s1';
mtj = 1;
figNamePrefix = 'none';
% figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\matlab_final\SOTA';
% figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\matlab_final\MLA_vs_WL';
% figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\matlab_final\best';

%%% select figures to make
makeplot.kinematics = 1;
makeplot.kinetics = 1;
makeplot.soleus = 0;
makeplot.sol_all = 0;
makeplot.GRF = 0;
makeplot.compareLiterature = 1;
makeplot.compareTakahashi17 = 0;
makeplot.compareTakahashi17_mtj_only = 1;
% makeplot.COP = 0;
makeplot.allQsTs = 0;
makeplot.k_mtj_lin = 0;
makeplot.windlass = 0;
makeplot.power = 0;
makeplot.power_T = 0;
makeplot.work = 0;
makeplot.work_bar = 0;
makeplot.power_main = 1;
makeplot.spatiotemp = 0;

%%% make separate fiures
if plot_full_separate
    PlotResults_3DSim_Report(ResultsFile,LegNames,RefData,mtj,makeplot,figNamePrefix);
end

%% make tabbed figure
% Set last argument to 0 for main figure, to 1 for cross-correlation
% coefficient figure or to 2 for both.
if plot_full_tabbed
    [h_default,h_ccc] = Plot3D(ResultsFile,reference_data,0);
end

%%
plot_COT_k = 0;
if plot_COT_k
%     idx = [1,5,9,13:16,21,24];
%     idx = [1,5,9];
%     idx = [1,13:16];   
idx=1;

    scs = get(0,'ScreenSize');
    figure('Position',[1+scs(3)/2,scs(4)/2+20,scs(3)/2, scs(4)/2-100]);
    CsV = hsv(length(idx)+1);
    load(ref{1},'R');
    subplot(1,3,1:2)
    plot([0,1e6],[1,1]*R.COT,'Color',CsV(1,:),'DisplayName','Falisse 2019')
    hold on
    grid on
    subplot(2,3,6)
    semilogx([10,2e6],[1,1]*R.COT,'Color',CsV(1,:),'DisplayName','Falisse 2019')
    hold on
    grid on
    ds = zeros(9,20);
    for i=1:length(idx)
        k_COT = zeros(4,numel(filteredResultsWithRef{idx(i)}));

        for j=1:numel(filteredResultsWithRef{idx(i)})
            load(filteredResultsWithRef{idx(i)}{j},'R');
            has_mtj = isfield(R.S,'mtj') && ~isempty(R.S.mtj) && R.S.mtj;
            if has_mtj
                k_COT(1,j) = R.S.kMT_li;
                k_COT(2,j) = R.COT;
                
                imtj = find(strcmp(R.colheaders.joints,'mtj_angle_r'));
                qdot_mtj = R.Qdots(:,imtj)*pi/180;
                M_li = R.windlass.M_li;
                P_mtj_li = qdot_mtj.*M_li/R.body_mass;
                P_mtj_li_pos = P_mtj_li;
                P_mtj_li_pos(P_mtj_li_pos<0) = 0;
                W_mtj_li = trapz(R.t,P_mtj_li_pos);
                k_COT(3,j) = max(P_mtj_li);
                k_COT(4,j) = W_mtj_li;
                
            end

        end
        subplot(1,3,1:2)
        plot(k_COT(1,2:end),k_COT(2,2:end),'.-','MarkerSize',20,'Color',CsV(1+i,:),'DisplayName',groupNames{idx(i)})
        subplot(2,3,6)
        semilogx(k_COT(1,2:end),k_COT(2,2:end),'.-','MarkerSize',10,'Color',CsV(1+i,:),'DisplayName',groupNames{idx(i)})
    end

    subplot(1,3,1:2)
    lh=legend('Location','northwest');
    lhPos = lh.Position;
    lhPos(1) = lhPos(1)+0.5;
    lhPos(2) = lhPos(2)+0.05;
    set(lh,'position',lhPos);
    xlabel('mtj stiffness (Nm/rad)')
    ylabel('COT (J kg^-^1 m^-^1)','Interpreter','tex')
    title('Cost of transport and midtarsal joint stiffness')
    xlim([0,2050])
    
    subplot(2,3,6)
    xlabel('\fontsize{10} mtj stiffness (Nm/rad)','Interpreter','tex')
    ylabel('\fontsize{10} COT (J kg^-^1 m^-^1)','Interpreter','tex')
    title('\fontsize{10} Convergence','Interpreter','tex')
    xlim([20,1e6])
    a1 = gca;
    a1.YAxisLocation = 'right';
    
end

%%
plot_frq_k = 0;
if plot_frq_k
%     idx = [1,5,9,13:16,19,21,24,25];
%     idx = [1,5,9];
%     idx = [1,13:16];   

    scs = get(0,'ScreenSize');
    figure('Position',[1+scs(3)/2,scs(4)/2+20,scs(3)/2, scs(4)/2-100]);
    CsV = hsv(length(idx)+1);
    load(ref{1},'R');
    subplot(1,3,1:2)
    plot([0,1e6],[1,1]*1./(R.tf_step*2),'Color',CsV(1,:),'DisplayName','Falisse 2019')
    hold on
    grid on
    subplot(2,3,6)
    semilogx([10,2e6],[1,1]*1./(R.tf_step*2),'Color',CsV(1,:),'DisplayName','Falisse 2019')
    hold on
    grid on
    ds = zeros(9,20);
    for i=1:length(idx)
        k_frq = zeros(2,numel(filteredResultsWithRef{idx(i)}));

        for j=1:numel(filteredResultsWithRef{idx(i)})
            load(filteredResultsWithRef{idx(i)}{j},'R');
            has_mtj = isfield(R.S,'mtj') && ~isempty(R.S.mtj) && R.S.mtj;
            if has_mtj
                k_frq(1,j) = R.S.kMT_li;
                k_frq(2,j) = 1./(R.tf_step*2);
            end

        end
        subplot(1,3,1:2)
        plot(k_frq(1,2:end),k_frq(2,2:end),'.-','MarkerSize',20,'Color',CsV(1+i,:),'DisplayName',groupNames{idx(i)})
        subplot(2,3,6)
        semilogx(k_frq(1,2:end),k_frq(2,2:end),'.-','MarkerSize',10,'Color',CsV(1+i,:),'DisplayName',groupNames{idx(i)})
    end

    subplot(1,3,1:2)
    lh=legend('Location','northwest');
    lhPos = lh.Position;
    lhPos(1) = lhPos(1)+0.5;
    lhPos(2) = lhPos(2)+0.05;
    set(lh,'position',lhPos);
    xlabel('mtj stiffness (Nm/rad)')
    ylabel('frequenct (Hz)','Interpreter','tex')
    title('Stride frequency and midtarsal joint stiffness')
    xlim([0,2050])
    
    subplot(2,3,6)
    xlabel('\fontsize{10} mtj stiffness (Nm/rad)','Interpreter','tex')
    ylabel('\fontsize{10} f (Hz)','Interpreter','tex')
    title('\fontsize{10} Convergence','Interpreter','tex')
    xlim([20,1e6])
    a1 = gca;
    a1.YAxisLocation = 'right';
    
end

%%
plot_aSol_k = 0;
if plot_aSol_k
%     idx = [1,5,9,13:16,19,21,24,25];
%     idx = [1,5,9];
%     idx = [1,13:16];   

    scs = get(0,'ScreenSize');
    figure('Position',[1+scs(3)/2,scs(4)/2+20,scs(3)/2, scs(4)/2-100]);
    CsV = hsv(length(idx)+1);
    load(ref{1},'R');
    iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
    subplot(1,3,1:2)
    plot([0,1e6],[1,1]*max(R.a(:,iSol)),'Color',CsV(1,:),'DisplayName','Falisse 2019')
    hold on
    grid on
    subplot(2,3,6)
    semilogx([10,2e6],[1,1]*max(R.a(:,iSol)),'Color',CsV(1,:),'DisplayName','Falisse 2019')
    hold on
    grid on
    ds = zeros(9,20);
    for i=1:length(idx)
        iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
        k_aSol = zeros(2,numel(filteredResultsWithRef{idx(i)}));

        for j=1:numel(filteredResultsWithRef{idx(i)})
            load(filteredResultsWithRef{idx(i)}{j},'R');
            has_mtj = isfield(R.S,'mtj') && ~isempty(R.S.mtj) && R.S.mtj;
            if has_mtj
                k_aSol(1,j) = R.S.kMT_li;
                k_aSol(2,j) = max(R.a(:,iSol));
            end

        end
        subplot(1,3,1:2)
        plot(k_aSol(1,2:end),k_aSol(2,2:end),'.-','MarkerSize',20,'Color',CsV(1+i,:),'DisplayName',groupNames{idx(i)})
        subplot(2,3,6)
        semilogx(k_aSol(1,2:end),k_aSol(2,2:end),'.-','MarkerSize',10,'Color',CsV(1+i,:),'DisplayName',groupNames{idx(i)})
    end

    subplot(1,3,1:2)
    lh=legend('Location','northwest');
    lhPos = lh.Position;
    lhPos(1) = lhPos(1)+0.5;
    lhPos(2) = lhPos(2)+0.05;
    set(lh,'position',lhPos);
    xlabel('mtj stiffness (Nm/rad)')
    ylabel('max activity','Interpreter','tex')
    title('Max soleus activity and midtarsal joint stiffness')
    xlim([0,2050])
    
    subplot(2,3,6)
    xlabel('\fontsize{10} mtj stiffness (Nm/rad)','Interpreter','tex')
    ylabel('\fontsize{10} a_max','Interpreter','tex')
    title('\fontsize{10} Convergence','Interpreter','tex')
    xlim([20,1e6])
    a1 = gca;
    a1.YAxisLocation = 'right';
    
end

%%
plot_COT_v = 0;
if plot_COT_v
    % data from http://dx.doi.org/10.1098/rsif.2019.0402
    x_Fal = 0.73:0.1:2.73;
    y_Fal = [4.8577,4.4555, 4.0769, 3.9148, 3.7499, 3.6131, 3.5464, 3.4393,...
        3.3582, 3.3896, 3.4203, 3.3646, 3.4244, 3.4749, 3.5139, 3.7447,...
        3.7519, 3.7449, 3.6630, 3.6883, 3.4600];

    % data from S. Song en H. Geyer, „The Energetic Cost of Adaptive Feet in Walking,” 
    % in International Conference on Robotics and Biomimetics, Phuket, Thailand, 2011. 
    x_Song = 0.8:0.2:1.8;
    y_Song_b = [4.55, 3.4, 2.95, 2.75, 2.85, 2.87];
    y_Song_h = [5.65, 4.12, 3.55, 2.9, 3, 2.73];

    for i=1:length(Results_vs)-1
        load(Results_vs{i},'R');
        x_Fal_mtp(i) = R.S.v_tgt;
        y_Fal_mtp(i) = R.COT;
    end
    for i=1:length(Results_vs_Song)-1
        load(Results_vs_Song{i},'R');
        x_Fal_mtj_Song(i) = R.S.v_tgt;
        y_Fal_mtj_Song(i) = R.COT;
    end

    scs = get(0,'ScreenSize');
    figure('Position',[1,scs(4)/2+20,scs(3)/2, scs(4)/2-100]);
    subplot(121)
    hold on
    grid on
    plot(x_Song,y_Song_b,'o-','DisplayName','Baseline (Song 2011)')
    plot(x_Song,y_Song_h,'d-','DisplayName','Human (Song 2011)')
    plot(x_Fal(1:11),y_Fal(1:11),'o-','DisplayName','Baseline (Falisse 2019)')
    plot(x_Fal_mtp,y_Fal_mtp,'v-','DisplayName','Baseline + mtp (Falisse 2019)')
    plot(x_Fal_mtj_Song,y_Fal_mtj_Song,'d-','DisplayName','Human (parameters Song 2011)')
    legend('Location','northeast')
    xlabel('v (m s^{-1})')
    ylabel('COT (J kg^{-1} m^{-1})')
    xlim([0.7,1.85])
    title('Cost of transport for different foot models')

    rel_diff_Song = (y_Song_b-y_Song_h)./y_Song_b*100;
    y_tmp = interp1(x_Fal,y_Fal,x_Fal_mtj_Song);
    rel_diff = (y_tmp-y_Fal_mtj_Song)./y_tmp*100;

    subplot(122)
    bar(x_Song,[rel_diff_Song;rel_diff]')
    hold on
    grid on
    legend({'Human (Song 2011) vs Baseline (Song 2011)','Human (parameters Song 2011) vs Baseline (Falisse 2019)'},'Location','southoutside')
    xlabel('v (m s^{-1})')
    ylabel('\delta COT (%)')
    xlim([0.7,1.85])
    title('Cost of transport relative to baseline model')
end

%%

if plot_static_foot
    
%     resultFiles = {fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Gefen2002_Gefen2002_Q-30_30_F0_3000_WLv3_ls150.mat'])
%                 fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_none_Gefen2002_Q-30_30_F0_3000_WLv3_ls150.mat'])
%                 fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_none_Ker1987_Q-30_30_F0_3000_WLv3_ls150.mat'])};
    
%     resultFiles = {fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_fitted4_Q-20_30_F0_0_WLv3_ls148_mtp1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_fitted4_Q-20_30_F0_0_WLv3_ls148_mtp2.mat'])};
    
%     resultFiles =  {fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k50_Q0_30_F0_1000_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k300_Q0_30_F0_1000_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k800_Q0_30_F0_1000_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_linear_k300_Q0_30_F0_1000_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Gefen2002_Gefen2002_Q0_30_F0_1000_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Song2011_Song2011_Q0_30_F0_1000_WLv3_ls150_sb1.mat'])};
    
%     resultFiles =  {fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k30_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k50_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k100_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k150_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k200_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k250_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k300_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k350_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k380_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k400_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k450_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k500_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k600_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k700_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k800_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k900_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k1000_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k1200_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k1500_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k2000_Q0_30_F0_960_WLv3_ls150_sb1.mat'])
%                     fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k5000_Q0_30_F0_960_WLv3_ls150_sb1.mat'])};


    resultFiles =  {fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k10_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k30_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k50_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k100_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k150_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k200_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k250_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k300_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k350_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k400_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k450_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k500_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k550_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k600_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k650_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k700_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k750_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k800_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k900_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k1000_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k1100_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k1200_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k1500_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k2000_Q0_30_F0_945_WLv3_ls150_sb1.mat'])
                    };
                
    % call plot function
    nrf = numel(resultFiles);
    CsV = hsv(nrf);
    for i=1:nrf
        load(resultFiles{i},'R');
        if i==1
            h = PlotResults_FootSim(R,CsV(i,:),0,6);
        else
            PlotResults_FootSim(R,CsV(i,:),h,6);
        end
    end
    
end















%%

% k_MT = [10,30,50:50:800,900:100:1200,1500,2000];
% mtp_qdkl = [-0.1419,-0.064616,-0.016848,0.04365,0.067051,0.075793,0.078264,0.077705,...
%     0.075768,0.073134,0.0703,0.067442,0.068885,0.061921,0.059326,0.056915,...
%     0.054708,0.052554,0.047983,0.045381,0.043247,0.046333,0.032625,0.026503];
% 
% scs = get(0,'ScreenSize');
% h1 = figure('Position',[scs(3)/2,140,scs(3)/2, scs(4)/2-100]);
% subplot(1,2,1)
% 
% plot(k_MT,mtp_qdkl*100,'.-','MarkerSize',20,'DisplayName','Relative mtp quasi-stiffness increase with load')
% ylabel('Relative mtp quasi-stiffness increase with load (%)')
% xlim([-50,2050])
%  
% yyaxis right
% plot(k_COT(1,2:end),k_COT(2,2:end),'.-','MarkerSize',20,'DisplayName','Cost of transport')
% ylabel('COT (J kg^-^1 m^-^1)')
% 
% set(gca,'YColor',[0,0,0])
% yyaxis left
% 
% mtp_qdkl_i = interp1(k_MT,mtp_qdkl,k_COT(1,2:12),'spline','extrap');
% rr1 = xcorr(mtp_qdkl_i,k_COT(2,2:12),0,'coeff');
% text(1500,-7,['R = ' num2str(rr1,2)])
% 
% legend('Location','northoutside')
% xlabel('Midtarsal stiffness')
% 
% 
% subplot(1,2,2)
% plot(k_MT,mtp_qdkl*100,'.-','MarkerSize',20,'DisplayName','Relative mtp quasi-stiffness increase with load')
% ylabel('Relative mtp quasi-stiffness increase with load (%)')
% xlim([-50,2050])
%  
% yyaxis right
% plot(k_COT(1,2:end),k_COT(4,2:end),'.-','MarkerSize',20,'DisplayName','Positive Work by k_{mtj}')
% ylabel('Pos W (J/kg)')
% 
% set(gca,'YColor',[0,0,0])
% yyaxis left
% 
% rr2 = xcorr(mtp_qdkl_i,k_COT(4,2:12),0,'coeff');
% text(1500,-7,['R = ' num2str(rr2,2)])
% 
% legend('Location','northoutside')
% xlabel('Midtarsal stiffness')
% 
% rr3 = xcorr(k_COT(2,2:12),k_COT(4,2:12),0,'coeff');
% 
% 
% figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\matlab_final\FootSim';
% set(h1,'PaperPositionMode','auto')
% print(h1,[figNamePrefix '_mtp_stiffening'],'-dpng','-r0')
% print(h1,[figNamePrefix '_mtp_stiffening'],'-depsc')



