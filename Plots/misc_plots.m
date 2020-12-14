

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



%% tmt stiffness

qd = linspace(-2,4,1000);
q=qd*pi/180;
k=800;
T1 = k*q;

figure
plot(qd,T1,'LineWidth',2)
hold on
grid on
xlabel('angle (°)')
ylabel('passive torque (Nm)')
title('Nonlinear tendon and ligament stiffness')

k1 = 800;
k2 = 2;
t1 = 0.5*pi/180;
T2 = k1.*(q-t1*tanh(q/(k2*t1)));
plot(qd,T2,'LineWidth',2)


% k1 = 800;
% k2 = 1;
% t1 = 1*pi/180;
% T2 = k1.*(q-t1*tanh(q/(k2*t1)));
% plot(qd,T2)
% 
% k1 = 800;
% k2 = 3;
% t1 = 0.5*pi/180;
% T2 = k1.*(q-t1*tanh(q/(k2*t1)));
% plot(qd,T2)
% 
% k1 = 800;
% k2 = 3;
% t1 = 1.5*pi/180;
% T2 = k1.*(q-t1*tanh(q/(k2*t1)));
% plot(qd,T2)


%%

figure
load('D:\school\WTK\thesis\model\3dpredictsim\Results\all\Pog_s1_bCst_act_pp.mat');
x = 1:(100-1)/(size(R.Qs,1)-1):100;
plot(x,R.Qs(:,14),'DisplayName','h = cst');
hold on
grid on
title('Knee angle (assistance)')
ylabel('Angle (°)')
xlabel('Gait cycle (%)')
load('D:\school\WTK\thesis\model\3dpredictsim\Results\all\Pog_s1_aCst_act_pp.mat');
plot(x,R.Qs(:,14),'--','DisplayName','\alpha = cst');
legend('location','best')

%%
figure
load('D:\school\WTK\thesis\model\3dpredictsim\Results\all\Pog_s1_bCst_pp.mat');
x = 1:(100-1)/(size(R.Qs,1)-1):100;
plot(x,R.Tid(:,16),'DisplayName','h = cst');
hold on
grid on
title('Ankle torque (normal shoes)')
ylabel('Torque (Nm)')
xlabel('Gait cycle (%)')
load('D:\school\WTK\thesis\model\3dpredictsim\Results\all\Pog_s1_aCst_pp.mat');
plot(x,R.Tid(:,16),'--','DisplayName','\alpha = cst');
legend('location','best')

%%
figure
load('D:\school\WTK\thesis\model\3dpredictsim\Results\all\Pog_s1_bCst_pas_pp.mat');
iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
x = 1:(100-1)/(size(R.Qs,1)-1):100;
plot(x,R.a(:,iSol),'DisplayName','h = cst');
hold on
grid on
title('Soleus activity (unpowered exo)')
ylabel('Activity (-)')
xlabel('Gait cycle (%)')
load('D:\school\WTK\thesis\model\3dpredictsim\Results\all\Pog_s1_aCst_pas_pp.mat');
plot(x,R.a(:,iSol),'--','DisplayName','\alpha = cst');
legend('location','best')
ylim([0.02, 0.36])



