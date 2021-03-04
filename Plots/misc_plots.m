

%% plot soleus activity Normal
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

%% plot soleus activity Active
load('D:\school\WTK\thesis\model\3dpredictsim\Data\Pog_s1.mat','Dat');
Qref = Dat.Active.gc;

iSol_data = find(strcmp(Dat.Active.EMGheaders,'soleus_r'));
sol_act_Active = [Dat.Active.gc.lowEMG_mean(ceil(end/2):end,iSol_data); Dat.Active.gc.lowEMG_mean(1:ceil(end/2)-1,iSol_data)];

idx_jref = strcmp(Qref.colheaders,'ankle_angle_r');
T_ankle = -Qref.Tall_bio_mean(:,idx_jref);
x1 = Dat.Active.Step.PercStance;

figure
% subplot(211)
plot(sol_act_Active)
hold on
% line([x1 x1], get(gca, 'ylim'),'color','k','LineStyle','--')
xlabel('Gait cycle (%)')
ylabel('Soleus activity')
% grid on
% subplot(212)
% plot(T_ankle)
% hold on
% % line([x1 x1], get(gca, 'ylim'),'color','k','LineStyle','--')
% xlabel('Gait cycle (%)')
% ylabel('Ankle torque (Nm)')
% % grid on

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

%% plot starting point for important joints

type = 'Normal';
Cs = 'b';
figure
sgtitle('Joint angles and torques (without exoskeleton)')
joints_ref = {'hip_flexion_r','hip_adduction_r','knee_angle_r','ankle_angle_r'};
        
 joints_tit = {'Pelvis tilt','Pelvis list','Pelvis rotation','Pelvis tx',...
        'Pelvis ty','Pelvis tz','Hip flexion L','Hip adduction L',...
        'Hip rotation L','Hip flexion R','Hip adduction R','Hip rotation R',...
        'Knee L','Knee R','Ankle L','Ankle R',...
        'Subtalar L','Subtalar R','TMT L','TMT R','MTP L','MTP R',...
        'Lumbar extension','Lumbar bending','Lumbar rotation',...
        'Arm flexion L','Arm adduction L','Arm rotation L',...
        'Arm flexion R','Arm adduction R','Arm rotation R',...
        'Elbow flexion L','Elbow flexion R'};

load('D:\school\WTK\thesis\model\3dpredictsim\Data\Pog_s1.mat','Dat');
Qref = Dat.(type).gc;
load('D:\school\WTK\thesis\model\3dpredictsim\Results\all\Pog_s1_bCst_pp.mat');

idx_Qs = [10,11,14,16];
idx_title = [10,11,14,16];

j = 0;
label_fontsize  = 12;
line_linewidth  = 2;
for i = 1:length(idx_title)
    subplot(2,4,i)
    x = 1:(100-1)/(size(R.Qs,1)-1):100;
    % Experimental data
    
    idx_jref = strcmp(Qref.colheaders,joints_ref{i});
    if sum(idx_jref) == 1
        meanPlusSTD = (Qref.Qall_mean(:,idx_jref) + 2*Qref.Qall_std(:,idx_jref)).*180/pi;
        meanMinusSTD = (Qref.Qall_mean(:,idx_jref) - 2*Qref.Qall_std(:,idx_jref)).*180/pi;
        stepQ = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
        intervalQ = 1:stepQ:size(R.Qs,1);
        sampleQ = 1:size(R.Qs,1);
        meanPlusSTD = interp1(intervalQ,meanPlusSTD,sampleQ);
        meanMinusSTD = interp1(intervalQ,meanMinusSTD,sampleQ);

        hold on
        fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k','DisplayName','MoCap');
        alpha(.25);
    end


    % Simulation results
    x = 1:(100-1)/(size(R.Qs,1)-1):100;
    hold on;
    j=j+1;
    if i == length(idx_title)
         plot(x,R.Qs(:,idx_Qs(j)),'color',Cs,'linewidth',line_linewidth,'DisplayName','sim');
    else
        plot(x,R.Qs(:,idx_Qs(j)),'color',Cs,'linewidth',line_linewidth);
    end

    % Plot settings

    set(gca,'Fontsize',label_fontsize);
    title(joints_tit{idx_title(i)},'Fontsize',label_fontsize);
    % Y-axis
    if i == 1
        ylabel('Angle (°)','Fontsize',label_fontsize);
    end
    % X-axis
    L = get(gca,'XLim');
    NumTicks = 2;

    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);

end

j = 0;
for i = 1:length(idx_title)
    subplot(2,4,i+4)
    x = 1:(100-1)/(size(R.Qs,1)-1):100;
    % Experimental data
    idx_jref = strcmp(Qref.colheaders,joints_ref{i});
    if sum(idx_jref) == 1
        meanPlusSTD = Qref.Tall_mean(:,idx_jref) + 2*Qref.Tall_std(:,idx_jref);
        meanMinusSTD = Qref.Tall_mean(:,idx_jref) - 2*Qref.Tall_std(:,idx_jref);
        stepID = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
        intervalID = 1:stepID:size(R.Qs,1);
        sampleID = 1:size(R.Qs,1);
        meanPlusSTD = interp1(intervalID,meanPlusSTD,sampleID);
        meanMinusSTD = interp1(intervalID,meanMinusSTD,sampleID);
        hold on
        fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k','DisplayName','MoCap');
        alpha(.25);
    end

    % Simulation results
    x = 1:(100-1)/(size(R.Qs,1)-1):100;
    hold on;
    j=j+1;
    if i == length(idx_title)
        plot(x,R.Tid(:,idx_Qs(j)),'color',Cs,'linewidth',line_linewidth,'DisplayName','Simulated');
    else
        plot(x,R.Tid(:,idx_Qs(j)),'color',Cs,'linewidth',line_linewidth);
    end
    % Plot settings

    % Plot settings
    set(gca,'Fontsize',label_fontsize);
    title(joints_tit{idx_title(i)},'Fontsize',label_fontsize);
    % Y-axis
    if i == 1
        ylabel('Torque (Nm)','Fontsize',label_fontsize);
    end
    % X-axis
    L = get(gca,'XLim');
    NumTicks = 2;
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        
end


lh=legend('-DynamicLegend','location','east');
lh.Interpreter = 'none';
lhPos = lh.Position;
%         lhPos(1) = lhPos(1)+0.2;
set(lh,'position',lhPos);

