clear
clc
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/VariousFunctions']);
AddCasadiPaths();
import casadi.*

%%
cd(['D:\school\WTK\thesis\model\3dpredictsim\CasADiFunctions\'...
    'casadi_s1Fal_MuscModel_bCst_PF_Gefen2002_ls150_MT_nl_Gefen2002_MTP_Mf1_PIM']);

f_lMT_vMT_dM_11 = Function.load('f_lMT_vMT_dM');
qin_11 = zeros(1,11);
[lMT_11,~,dM_11] = f_lMT_vMT_dM_11(qin_11,qin_11);


%%
cd(['D:\school\WTK\thesis\model\3dpredictsim\CasADiFunctions\'...
    'casadi_s1Fal_MuscModel_bCst_PF_Natali2010_ls150_MT_k300_MTP_T1_PFx10']);

f_lMT_vMT_dM_10 = Function.load('f_lMT_vMT_dM');
qin_10 = zeros(1,10);
[lMT_10,~,dM_10] = f_lMT_vMT_dM_10(qin_10,qin_10);

cd(pathHere);
%%
muscleNames = {'glut_med1_r','glut_med2_r','glut_med3_r',...
    'glut_min1_r','glut_min2_r','glut_min3_r','semimem_r',...
    'semiten_r','bifemlh_r','bifemsh_r','sar_r','add_long_r',...
    'add_brev_r','add_mag1_r','add_mag2_r','add_mag3_r','tfl_r',...
    'pect_r','grac_r','glut_max1_r','glut_max2_r','glut_max3_r',......
    'iliacus_r','psoas_r','quad_fem_r','gem_r','peri_r',...
    'rect_fem_r','vas_med_r','vas_int_r','vas_lat_r','med_gas_r',...
    'lat_gas_r','soleus_r','tib_post_r','flex_dig_r','flex_hal_r',...
    'tib_ant_r','per_brev_r','per_long_r','per_tert_r','ext_dig_r',...
    'ext_hal_r','ercspn_r','intobl_r','extobl_r','ercspn_l',...
    'intobl_l','extobl_l'};


%%
diff = full(lMT_10) - full(lMT_11);
diff2 = diff;
diff2(abs(diff)<=3e-3) = 0;

for i=1:length(diff)
    if diff2(i)
        disp([muscleNames{i} '   ' num2str(diff2(i)) '   ' num2str(full(lMT_10(i)/lMT_11(i)))])
    end
end
disp('  ')
disp('  ')

% [[dM_10(1:6);0;dM_10(7:end)],dM_11]

%%
ExtPoly = '_mtp';
subject = 'subject1';
pathmusclemodel = fullfile(pathRepo,'MuscleModel',subject);
load([pathmusclemodel,'/MTparameters_',subject, ExtPoly, '.mat']);
musi = 1:length(muscleNames(1:end-3));

params = MTparameters;

lMo = params(2,:);
lTs = params(3,:);
alphao = params(4,:);

lMTo = lMo.*cos(alphao) + lTs;

%%

for i=1:length(diff)
    if diff2(i)
        MTparameters(3,i) = MTparameters(3,i) - diff2(i);
        disp([muscleNames{i} '  ' num2str( lTs(i) ) '  ' num2str( MTparameters(3,i) ) ])
    end
end
disp('  ')
disp('  ')


lMo2 = MTparameters(2,:);
lTs2 = MTparameters(3,:);
alphao2 = MTparameters(4,:);

lMTo2 = lMo2.*cos(alphao2) + lTs2;

for i=1:length(diff)
    if diff2(i)
        rel10 = full(lMT_10(i))/lMTo(i);
        rel11 = full(lMT_11(i))/lMTo2(i);
        disp([muscleNames{i} '  ' num2str(rel10) '  ' num2str(rel11) '  ' num2str(rel11/rel10) ])
        
    end
end

%%

% ExtPoly = '_mtj';
% pathmusclemodel = fullfile(pathRepo,'MuscleModel',subject);
% save([pathmusclemodel,'/MTparameters_',subject, ExtPoly, '.mat'],'MTparameters');




















