


pathmain        = pwd;
[pathRepo,~,~]  = fileparts(pathmain);
folder = '\Figures';
file = 'arch_stiffness_Ker87.png';
pathRefImg = fullfile(pathRepo,folder,file);
img_Ker = imread(pathRefImg);
   
CsV = hsv(3);

figure
hold on
hi1 = image([1,9.65]*0.9,flip([0,4]*0.9^2),img_Ker);
uistack(hi1,'bottom')


Result_compl_st = load(fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Gefen2002_Gefen2002_Q-30_30_F0_3000_WLv3_ls150.mat']),'R');
R = Result_compl_st.R;
Fs_tib = R.Fs_tib;
l_fa = R.l_fa((R.Qs_mtp(:)==0),:);
plot((l_fa-l_fa(1))*1000,Fs_tib/1000,'-o','Color',CsV(3,:),'DisplayName','Compliant, intact foot')



Result_compl_noPF_st = load(fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_none_Gefen2002_Q-30_30_F0_3000_WLv3_ls150.mat']),'R');
R = Result_compl_noPF_st.R;
Fs_tib = R.Fs_tib;
l_fa = R.l_fa((R.Qs_mtp(:)==0),:);
plot((l_fa-R.L0)*1000,Fs_tib/1000,'--d','Color',CsV(3,:),'DisplayName','Compliant, plantar fascia removed')



Result_stiff_st = load(fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k300_Q-30_30_F0_3000_WLv3_ls150_sb1_PFx10.mat']),'R');
R = Result_stiff_st.R;
Fs_tib = R.Fs_tib;
l_fa = R.l_fa((R.Qs_mtp(:)==0),:);
plot((l_fa-l_fa(1))*1000,Fs_tib/1000,'-o','Color',CsV(2,:),'DisplayName','Combination, intact foot')


Result_stiff_noPF_st = load(fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_none_k300_Q-30_30_F0_3000_WLv3_ls150_sb1_PFx10.mat']),'R');
R = Result_stiff_noPF_st.R;
Fs_tib = R.Fs_tib;
l_fa = R.l_fa((R.Qs_mtp(:)==0),:);
plot((l_fa-R.L0)*1000,Fs_tib/1000,'--d','Color',CsV(2,:),'DisplayName','Combination, plantar fascia removed')

xlabel('Horizontal elongation (mm)')
ylabel('Vertical force (kN)')
title({'Foot arch stiffness','\rm as defined by Ker and all, 1987'})
lg12 = legend('Location','northwest');



