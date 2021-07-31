close all
clear
clc

load('D:\school\WTK\thesis\model\3dpredictsim\Results\Final\Fal_s1_bCst_ig1_pp.mat','R');
load('D:\school\WTK\thesis\model\3dpredictsim\Data\Fal_s1.mat','Dat');

iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
iGas = find(strcmp(R.colheaders.muscles,'lat_gas_r'));
iGas2 = find(strcmp(R.colheaders.muscles,'med_gas_r'));
if isempty(iGas)
    iGas = find(strcmp(R.colheaders.muscles,'gaslat_r'));
end
if isempty(iGas2)
    iGas2 = find(strcmp(R.colheaders.muscles,'gasmed_r'));
end

h3 = figure;

imus = [iSol,iGas,iGas2];

CsV = [0.3;0.5;0]';
inr = 1;
label_fontsize  = 12;
line_linewidth  = 0.5;
NumTicks = 6;
LegName = 'Fal';


%%
type = 'Normal';

iSol_data = find(strcmp(Dat.(type).EMGheaders,'Soleus'));
iGas_data = find(strcmp(Dat.(type).EMGheaders,'Gastrocnemius-medialis'));
iGas2_data = find(strcmp(Dat.(type).EMGheaders,'Gastrocnemius-lateralis'));

ankle_act(:,1) = Dat.(type).gc.lowEMG_mean([51:end,1:50],iSol_data);
ankle_act(:,2) = Dat.(type).gc.lowEMG_mean([51:end,1:50],iGas_data);
ankle_act(:,3) = Dat.(type).gc.lowEMG_mean([51:end,1:50],iGas2_data);

ankle_a = [ankle_act(ceil(end/2):end,:); ankle_act(1:ceil(end/2)-1,:)];

for imu=1:3
    subplot(5,3,imu); hold on;
    yyaxis right
    plot(ankle_a(:,imu),'-k','DisplayName','EMG data')
    a1 = gca;
    a1.YColor = [0,0,0];
    if imu==3
        ylabel('EMG (mV)')
    end
    yyaxis left
    a1 = gca;
    a1.YColor = [0,0,0];
end

%%

subplot(5,3,1); hold on;
plot(R.a(:,iSol),'-','Color',CsV(inr,:),'DisplayName',LegName);
title('Soleus')
ylabel('Activity (-)')
grid on
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

subplot(5,3,2); hold on;
plot(R.a(:,iGas),'-','Color',CsV(inr,:),'DisplayName',LegName);
title('Gastrocnemius-medialis')
grid on
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

subplot(5,3,3); hold on;
plot(R.a(:,iGas2),'-','Color',CsV(inr,:),'DisplayName',LegName);
title('Gastrocnemius-lateralis')
grid on
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))

% lh=legend('-DynamicLegend','location','northeast');
% lh.Interpreter = 'none';

for imu=1:3

    subplot(5,3,3+imu)
    plot(R.lMT(:,imus(imu)),'-','Color',CsV(inr,:)); hold on;
    hold on
    grid on
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    if imu==1
        ylabel('Length (m)');
    end

    subplot(5,3,6+imu)
    plot(R.lT(:,imus(imu)),'-','Color',CsV(inr,:)); hold on;
    hold on
    grid on
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    if imu==1
        ylabel('Tendon length (m)');
    end

    subplot(5,3,9+imu)
    plot(R.lMtilde(:,imus(imu)),'-','Color',CsV(inr,:)); hold on;
    hold on
    grid on
    L = get(gca,'XLim');
    set(gca,'XTick',linspace(L(1),L(2),NumTicks))
    if imu==1
        ylabel('Fibre length (-)');
    end
    ylim([0.6,1.05]);
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    
end
% lhPos = lh.Position;
% lhPos(1) = lhPos(1)+0.1;
% % lhPos(2) = lhPos(2)+0.08;
% set(lh,'position',lhPos);

iankle = strcmp(R.colheaders.joints,'ankle_angle_r');
iknee = strcmp(R.colheaders.joints,'knee_angle_r');

subplot(5,3,13)
plot(R.Qs(:,iankle),'linewidth',line_linewidth,'Color',CsV(inr,:));
title('ankle');
grid on
% axis tight
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
ylabel('Angle (°)');

subplot(5,3,14)
plot(R.Qs(:,iknee),'linewidth',line_linewidth,'Color',CsV(inr,:));
title('knee');
grid on
% axis tight
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
ylim([-60,10]);

% subplot(5,3,15)
% plot(R.Qs(:,iknee)+R.Qs(:,iankle),'linewidth',line_linewidth,'Color',CsV(inr,:));
% title('knee + ankle')
% grid on
% axis tight
% set(gca,'XTick',linspace(L(1),L(2),NumTicks))

