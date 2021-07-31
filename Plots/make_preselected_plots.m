clear all
close all
clc

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
ref{1} = fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_pp.mat']);
reference_data = 'norm';

%%


Results_default = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_MTP_T5_ig24_v2_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k150_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_MTP_T5_ig24_pp.mat'])
%     fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_MTP_T5_ig21_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k1500_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k2000_MTP_T5_ig24_pp.mat'])};

Results_d010 = {
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

Results_d020 = {
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

Results_d050 = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_d050_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_d050_MTP_T5_ig24_pp.mat'])
%     fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k250_d050_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_d050_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_d050_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_d050_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_d050_MTP_T5_ig24_pp.mat'])};

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

Results_vs = {
    fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_v8_pp.mat'])
    fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_v10_pp.mat'])
    fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_v12_pp.mat'])
    fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_v14_pp.mat'])
    fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_v16_pp.mat'])
    fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_v18_pp.mat'])
    fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig1_v27_pp.mat'])};
    
Results_vs_Song = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu0_ig1_v8_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu0_ig1_v10_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu0_ig1_v12_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu0_ig1_v14_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu0_ig1_v16_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu0_ig1_v18_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu0_ig1_v27_pp.mat'])};


Results_PFonly = {
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k10_MTP_T5_ig24_PFx10_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k10_MTP_T2_ig24_PFx10_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_T5_ig24_PFx10_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_T2_ig24_PFx10_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k50_MTP_T5_ig24_PFx5_pp.mat'};

Results_Musc = {
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k20_MTP_Mu5_ig24_PFx5_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu1_ig23_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu1_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu1_np_ig23_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu1_spx2_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu5_ig23_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu5_np_ig23_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_Mu2_ig24_PFx10_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_Mu5_ig24_PFx10_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k30_MTP_Mu5_ig24_PFx5_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_fitted6_MTP_Mu5_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_Mu5_ig24_pp.mat'};
    
Results_singed_lin = {
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_k100_10_MTP_T5_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_k200_10_MTP_T5_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_k300_10_MTP_T5_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_k400_10_MTP_T5_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_k500_10_MTP_T5_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_k800_10_MTP_T5_ig24_pp.mat'};

Results_singed_lin_test = {
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_T5_ig1_300_50_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_T5_ig24_300_10_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_T5_ig24_300_30_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_T5_ig24_500_30_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_T5_ig24_500_50_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_aCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_T5_ig24_500_50_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_Mu5_ig24_300_10_pp.mat'};

Results_other_PF = {
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Gefen2002_ls150_MT_nl_signed_lin_MTP_T5_ig24_300_10_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_Song2011_MTP_T5_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_signed_lin_MTP_T5_ig24_300_10_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Song2011_ls150_MT_nl_signed_lin_MTP_T5_ig24_800_50_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_linear_ls150_MT_nl_signed_lin_MTP_T5_ig24_300_10_pp.mat'};

Results_no_PF = {
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_None_ls150_MT_k250_MTP_T17_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_None_ls150_MT_k500_MTP_T17_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_None_ls150_MT_k800_MTP_T17_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_None_ls150_MT_k1000_MTP_T17_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_None_ls150_MT_k1500_MTP_T17_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_None_ls150_MT_k2000_MTP_T17_ig24_pp.mat'};
    
Results_mtp_d0 = {
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k100_MTP_T5_d00_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k200_MTP_T5_d00_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_d00_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k400_MTP_T5_d00_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k500_MTP_T5_d00_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k800_MTP_T5_d00_ig24_pp.mat'};

Results_stiff_contact = {
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_ig24_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_spx10_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu1_np_ig23_pp.mat'])
    fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_Mu1_spx2_ig24_pp.mat'])};

Results_misc = {
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_fitted6_MTP_T2_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_fitted6_MTP_T5_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_none_ls148_MT_k1000_MTP_k10_ig24_pp.mat'
    'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_none_ls148_MT_k500_MTP_k10_ig24_pp.mat'
    };



%% check for unused results
% Results_all = {Results_default{:}, Results_d010{:}, Results_d020{:}, Results_d050{:}, Results_T2{:},...
%     Results_T10{:}, Results_PFx2{:}, Results_vs_Song{:}, Results_PFonly{:}, Results_Musc{:},...
%     Results_singed_lin{:}, Results_singed_lin_test{:}, Results_other_PF{:}, Results_no_PF{:},...
%     Results_mtp_d0{:}, Results_stiff_contact{:}, Results_misc{:}};
% 
% dpath = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint';
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

%%


filteredResultsWithRef{1} = {ref{:}, Results_default{:}};               % 1
filteredResultsWithRef{end+1} = {ref{:}, Results_d010{:}};              % 2
filteredResultsWithRef{end+1} = {ref{:}, Results_d020{:}};              % 3
filteredResultsWithRef{end+1} = {ref{:}, Results_d050{:}};              % 4
filteredResultsWithRef{end+1} = {ref{:}, Results_T2{:}};                % 5
filteredResultsWithRef{end+1} = {ref{:}, Results_T10{:}};               % 6
filteredResultsWithRef{end+1} = {ref{:}, Results_PFx2{:}};              % 7
filteredResultsWithRef{end+1} = {ref{:}, Results_vs{:}};                % 8
filteredResultsWithRef{end+1} = {ref{:}, Results_vs_Song{:}};           % 9
filteredResultsWithRef{end+1} = {ref{:}, Results_PFonly{:}};            % 10
filteredResultsWithRef{end+1} = {ref{:}, Results_Musc{:}};              % 11
filteredResultsWithRef{end+1} = {ref{:}, Results_singed_lin{:}};        % 12
filteredResultsWithRef{end+1} = {ref{:}, Results_singed_lin_test{:}};   % 13
filteredResultsWithRef{end+1} = {ref{:}, Results_other_PF{:}};          % 14
filteredResultsWithRef{end+1} = {ref{:}, Results_no_PF{:}};             % 15
filteredResultsWithRef{end+1} = {ref{:}, Results_mtp_d0{:}};            % 16
filteredResultsWithRef{end+1} = {ref{:}, Results_stiff_contact{:}};     % 17
filteredResultsWithRef{end+1} = {ref{:}, Results_misc{:}};              % 18

groupNames = {'d_m_t_j = 0 Nms/rad, k_m_t_p = 5 Nm/rad','d_m_t_j = 1 Nms/rad','d_m_t_j = 2 Nms/rad','d_m_t_j = 5 Nms/rad',...
        'k_m_t_p = 2 Nm/rad','k_m_t_p = 10 Nm/rad','k_P_F x2','vs Falisse 2019','vs Song 2011','PF only','Muscle-driven mtp',...
        'k_m_t_j_,_n_e_g = 10 Nm/rad','k_m_t_j_,_n_e_g test','other PF models','no PF','d_m_t_p = 0 Nms/rad','stiffer contact','misc'};


%%
% [h_default,~] = Plot3D(filteredResultsWithRef{4},reference_data,0);

%%
ResultsFile = filteredResultsWithRef{1};
LegNames = {'Falisse 2019'};

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
makeplot.power = 0;
makeplot.power_T = 0;
makeplot.work = 1;

figNamePrefix = 'none';

% PlotResults_3DSim_Report(ResultsFile,LegNames,RefData,mtj,makeplot,figNamePrefix);


%%
plot_COT_k = 1;
if plot_COT_k
    scs = get(0,'ScreenSize');
    figure('Position',[1+scs(3)/2,scs(4)/2+20,scs(3)/2, scs(4)/2-100]);
    idx = [1:7,12,16];
    CsV = hsv(length(idx)+1);
    load(ref{1},'R');
    plot([0,2000],[1,1]*R.COT,'Color',CsV(1,:),'DisplayName','Falisse 2019')
    hold on
    grid on
    ds = zeros(9,20);
    for i=1:length(idx)
        k_COT = zeros(2,numel(filteredResultsWithRef{idx(i)}));

        for j=1:numel(filteredResultsWithRef{idx(i)})
            load(filteredResultsWithRef{idx(i)}{j},'R');
            has_mtj = isfield(R.S,'mtj') && ~isempty(R.S.mtj) && R.S.mtj;
            if has_mtj
                k_COT(1,j) = R.S.kMT_li;
                k_COT(2,j) = R.COT;
            end

        end
        plot(k_COT(1,2:end),k_COT(2,2:end),'.-','Color',CsV(1+i,:),'DisplayName',groupNames{idx(i)})
    end

    legend('Location','northeastoutside')
    xlabel('mtj stiffness (Nm/rad)')
    ylabel('COT (J/(kg m)')
    title('Cost of transport and midtarsal joint stiffness')
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

