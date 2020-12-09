

%% plot soleus activity
load('D:\school\WTK\thesis\model\3dpredictsim\Data\Pog_s1.mat','Dat');
Qref = Dat.Normal.gc;

iSol_data = find(strcmp(Dat.Normal.EMGheaders,'soleus_r'));
sol_act_normal = [Dat.Normal.gc.lowEMG_mean(ceil(end/2):end,iSol_data); Dat.Normal.gc.lowEMG_mean(1:ceil(end/2)-1,iSol_data)];

idx_jref = strcmp(Qref.colheaders,'ankle_angle_r');
T_ankle = -Qref.Tall_bio_mean(:,idx_jref);
x1 = Dat.Normal.Step.PercStance;

figure
subplot(211)
plot(sol_act_normal)
hold on
% line([x1 x1], get(gca, 'ylim'),'color','k','LineStyle','--')
xlabel('Gait cycle (%)')
ylabel('Soleus activity')
grid on
subplot(212)
plot(T_ankle)
hold on
% line([x1 x1], get(gca, 'ylim'),'color','k','LineStyle','--')
xlabel('Gait cycle (%)')
ylabel('Ankle torque (Nm)')
grid on

%% plot assistance profile

load('D:\school\WTK\thesis\model\3dpredictsim\Results\Batchsim_tmt_linear\Pog_s1_tmt_bCst_d00_k800_ig24_act_pp.mat');
figure
plot(-R.T_exo(:,2))
ylabel('Exo Moment - Ankle [Nm]')
xlabel('% Gait cycle')
title('Assistance profile')



