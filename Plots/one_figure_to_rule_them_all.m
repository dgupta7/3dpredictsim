clear all
close all
clc

[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);

%%
col_1 = [231,176,55]/256;
col_2 = [212,216,66]/256;
col_3 = [82,189,236]/256;
col_4 = [29,141,176]/256;
col_5 = [47,77,93]/256;

cols = {col_2,col_1,col_3,col_4,col_5};
% cols = {col_2,col_3};

scs = get(0,'ScreenSize');
f1 = figure('Position',[10-scs(3),50,scs(3)/2, (scs(4)-140)]);
nh = 4;
nv = 5;
label_fontsize = 10;

set(f1,'Color','w');
%%
subplot(nv,nh,3)
lh=legend('location','northwest');
lhPos = lh.Position;
lhPos(1) = lhPos(1)+0.01;
lhPos(2) = lhPos(2)-0.2;
set(lh,'position',lhPos);
% title(lh,'Dynamische simulatie')
lh.Box = 'off';
LegName = {'Voet met 2 segmenten','Voet met 3 segmenten (soepel)','Voet met 3 segmenten (stug)'};



%%
% subplot(nv,nh,7)
% hold on
% axis equal
% plot([27,0,5,7],[0,0,10,20],'Color',cols{1},'LineWidth',1)
% plot([20,5],[0,10],'.-','Color',cols{1},'LineWidth',1,'MarkerSize',15)
% xlim([-3,30])
% tmp = gca;
% tmp.YAxis.Visible = 'off';
% tmp.XAxis.Visible = 'off';
% title('\rm 2 segmenten')
% 
% subplot(nv,nh,8)
% hold on
% axis equal
% plot([27,20,10,0,5,7],[0,0,6,0,10,20],'Color',cols{2},'LineWidth',1)
% plot([20,10,5],[0,6,10],'.-','Color',cols{2},'LineWidth',1,'MarkerSize',15)
% xlim([-3,30])
% tmp = gca;
% tmp.YAxis.Visible = 'off';
% tmp.XAxis.Visible = 'off';
% title('\rm 3 segmenten')

%%
subplot(nv,nh,[15:16]);
pathRefImg = fullfile(pathRepo,'\Figures\gait_cycle.png');
img_GC = imread(pathRefImg);
hold on
% axis tight
axis equal
hi3 = image([-5,65],flip([0,12.75]),img_GC);
uistack(hi3,'bottom')
xlim([-5,100])
ylim([-1,14])
xlabel('Stapcyclus (%)')
ax1=gca;
ax1.XTick = [0:10:100];
ax1.YTick = '';
ax1.YAxis.Visible = 'off';

% ax2 = axes('Position', get(ax1,'Position'),'XAxisLocation','top','Color','none','XColor','k');
% ax2.XTick = [0:10:100];
% ax2.YTick = '';
% ax2.XLim = ax1.XLim/0.6;
% axis equal
% xlabel('Stand (%)')
% ax2.YAxis.Visible = 'off';


%%




%%

Result_baseline = load(fullfile([pathRepo '\Results\Final\Fal_s1_bCst_ig21_pp.mat']),'R');
Result_stiff = load(fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T1_spx10_ig23_PFx10_spx10_pp.mat']),'R');
Result_compl = load(fullfile([pathRepo '\Results\MidTarsalJoint\Fal_s1_bCst_PF_Gefen2002_ls150_MT_nl_Gefen2002_MTP_T5_ig23_pp.mat']),'R');
Results_PredSim = {Result_baseline,Result_compl,Result_stiff};

%%
R = Results_PredSim{2}.R;
x = 1:(100-1)/(size(R.Qs,1)-1):100;

load('D:\school\WTK\thesis\model\3dpredictsim\Data\Fal_s1.mat','Dat');
Qref = Dat.Normal.gc;

ihip1_ref = strcmp(Qref.colheaders,'hip_flexion');
% ihip2_ref = strcmp(Qref.colheaders,'hip_adduction');
iknee_ref = strcmp(Qref.colheaders,'knee_angle');
iankle_ref = strcmp(Qref.colheaders,'ankle_angle');
i_ref = [ihip1_ref;iknee_ref;iankle_ref];

for i=1:3
    subplot(nv,nh,i)
    hold on
    idx_jref = i_ref(i,:);
    meanPlusSTD = (Qref.Qall_mean(:,idx_jref) + 2*Qref.Qall_std(:,idx_jref)).*180/pi;
    meanMinusSTD = (Qref.Qall_mean(:,idx_jref) - 2*Qref.Qall_std(:,idx_jref)).*180/pi;
    stepQ = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
    intervalQ = 1:stepQ:size(R.Qs,1);
    sampleQ = 1:size(R.Qs,1);
    meanPlusSTD = interp1(intervalQ,meanPlusSTD,sampleQ);
    meanMinusSTD = interp1(intervalQ,meanMinusSTD,sampleQ); 
    fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],[0.7,0.7,0.7],'DisplayName','Referentiedata (\pm2 std) [1]');
    alpha(.25);
    set(gca,'XTick',[0:20:100])
    xlabel('Stapcyclus (%)','Fontsize',label_fontsize);
end
subplot(nv,nh,1)
ylabel('Positie (°)','Fontsize',label_fontsize);
subplot(nv,nh,nh+1)
ylabel('Positie (°)','Fontsize',label_fontsize);

for i=1:3
    subplot(nv,nh,2*nh+i)
    hold on
    idx_jref = i_ref(i,:);
    meanPlusSTD = Qref.Tall_mean(:,idx_jref) + 2*Qref.Tall_std(:,idx_jref);
    meanMinusSTD = Qref.Tall_mean(:,idx_jref) - 2*Qref.Tall_std(:,idx_jref);
    stepQ = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
    intervalQ = 1:stepQ:size(R.Qs,1);
    sampleQ = 1:size(R.Qs,1);
    meanPlusSTD = interp1(intervalQ,meanPlusSTD,sampleQ);
    meanMinusSTD = interp1(intervalQ,meanMinusSTD,sampleQ);
    fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],[0.7,0.7,0.7],'DisplayName','Experimentele data (\pm 2 std) [1]');
    alpha(.25);
    xlabel('Stapcyclus (%)','Fontsize',label_fontsize);
    set(gca,'XTick',[0:20:100])
end
subplot(nv,nh,2*nh+1)
ylabel('Moment (Nm)','Fontsize',label_fontsize);
subplot(nv,nh,3*nh+1)
ylabel('Moment (Nm)','Fontsize',label_fontsize);


%%
% subplot(nv,nh,2*nh)
% iSol_data = find(strcmp(Dat.Normal.EMGheaders,'Soleus'));
% 
% hold on
% yyaxis right
% meanPlusSTD = Dat.Normal.gc.lowEMG_mean(:,iSol_data) + 2*Dat.Normal.gc.lowEMG_std(:,iSol_data);
% meanMinusSTD = Dat.Normal.gc.lowEMG_mean(:,iSol_data) - 2*Dat.Normal.gc.lowEMG_std(:,iSol_data);
% stepQ = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
% intervalQ = 1:stepQ:size(R.Qs,1);
% sampleQ = 1:size(R.Qs,1);
% meanPlusSTD = interp1(intervalQ,meanPlusSTD,sampleQ);
% meanMinusSTD = interp1(intervalQ,meanMinusSTD,sampleQ);
% fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k','DisplayName','Experimentele data (\pm 2 std) [1]');
% alpha(.25);
% set(gca,'XTick',[0:20:100])
%     
% % plot(Dat.Normal.gc.lowEMG_mean(:,iSol_data),'-k','DisplayName','EMG data')
% a1 = gca;
% a1.YColor = [0,0,0];
% ylabel('EMG (mV)')
% axis tight
% yl = get(gca, 'ylim');
% ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)])
% xlim([0,100])
% yyaxis left
% a1 = gca;
% a1.YColor = [0,0,0];


%%

subplot(nv,nh,3*nh)
hold on
meanPlusSTD = Dat.Normal.gc.GRF.Fmean(:,2) + 2*Dat.Normal.gc.GRF.Fstd(:,2);
meanMinusSTD = Dat.Normal.gc.GRF.Fmean(:,2) - 2*Dat.Normal.gc.GRF.Fstd(:,2);
stepQ = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
intervalQ = 1:stepQ:size(R.Qs,1);
sampleQ = 1:size(R.Qs,1);
meanPlusSTD = interp1(intervalQ,meanPlusSTD,sampleQ);
meanMinusSTD = interp1(intervalQ,meanMinusSTD,sampleQ);
fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],[0.7,0.7,0.7],'DisplayName','Experimentele data (\pm 2 std) [1]');
alpha(.25);
set(gca,'XTick',[0:20:100])
set(gca,'YAxisLocation','right')


%%
for i=1:length(Results_PredSim)
    R = Results_PredSim{i}.R;
    ihip1 = find(strcmp(R.colheaders.joints,'hip_flexion_r'));
    ihip2 = find(strcmp(R.colheaders.joints,'hip_adduction_r'));
    iknee = strcmp(R.colheaders.joints,'knee_angle_r');
    iankle = strcmp(R.colheaders.joints,'ankle_angle_r');
    isubt = strcmp(R.colheaders.joints,'subtalar_angle_r');
    imtj = find(strcmp(R.colheaders.joints,'mtj_angle_r'));
    imtp = find(strcmp(R.colheaders.joints,'mtp_angle_r'));
    iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));


    %%
    P_ankle = R.Qdots(:,iankle)*pi/180.*R.Tid(:,iankle)/R.body_mass;
    P_subt = R.Qdots(:,isubt)*pi/180.*R.Tid(:,isubt)/R.body_mass;
    P_mtp = R.Qdots(:,imtp)*pi/180.*R.Tid(:,imtp)/R.body_mass;
    if ~isempty(imtj)
        P_mtj = R.Qdots(:,imtj)*pi/180.*R.Tid(:,imtj)/R.body_mass;
    else
        P_mtj = zeros(size(P_mtp));
    end
    P_joints = P_mtj + P_ankle + P_subt + P_mtp;
    
%     P_HC_heel = R.P_mech_contact.vertical.calcn.r/R.body_mass;
%     P_HC_ball = R.P_mech_contact.vertical.metatarsi.r/R.body_mass;
%     P_HC_toes = R.P_mech_contact.vertical.toes.r/R.body_mass;
%     P_HC = P_HC_heel + P_HC_ball + P_HC_toes;
%     P_tot = P_HC + P_joints;
    
    
    %%
    subplot(nv,nh,1)
    title('\rm Heup')
    hold on
    plot(x,R.Qs(:,ihip1),'Color',cols{i})
    axis tight
    yl = get(gca, 'ylim');
    ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
    xlim([0,100])
    set(gca,'XTick',[0:20:100])
    
%     subplot(nv,nh,2)
%     title('\rm Heup (adductie)')
%     hold on
%     plot(x,R.Qs(:,ihip2),'Color',cols{i})
%     axis tight
%     xlim([0,100])
        
    subplot(nv,nh,2)
    title({'Dynamische simulatie','\rm ','\rm Knie'})
    hold on
    plot(x,R.Qs(:,iknee),'Color',cols{i},'LineWidth',1)
    axis tight
    yl = get(gca, 'ylim');
    ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
    xlim([0,100])
    set(gca,'XTick',[0:20:100])
%     set(gca,'YTick',[-60,-30,0])
    
    subplot(nv,nh,3)
    title('\rm Enkel')
    hold on
    plot(x,R.Qs(:,iankle),'Color',cols{i},'LineWidth',1,'DisplayName',LegName{i})
    axis tight
    yl = get(gca, 'ylim');
    ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
    xlim([0,100])
    set(gca,'XTick',[0:20:100])
    
    if ~isempty(imtj)
        subplot(nv,nh,5)
        title('\rm Midtarsaal gewricht')
        hold on
        plot(x,R.Qs(:,imtj),'Color',cols{i},'LineWidth',1)
        axis tight
        yl = get(gca, 'ylim');
        ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
        xlim([0,100])
        set(gca,'XTick',[0:20:100])
        xlabel('Stapcyclus (%)','Fontsize',label_fontsize);
    end
    
    subplot(nv,nh,6)
    title('\rm Teenbasisgewricht')
    hold on
    plot(x,R.Qs(:,imtp),'Color',cols{i},'LineWidth',1)
    axis tight
    yl = get(gca, 'ylim');
    ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
    xlim([0,100])
    set(gca,'XTick',[0:20:100])
    xlabel('Stapcyclus (%)','Fontsize',label_fontsize);
    
    %%
    subplot(nv,nh,2*nh+1)
    title('\rm Heup')
    hold on
    plot(x,R.Tid(:,ihip1),'Color',cols{i},'LineWidth',1)
    axis tight
    yl = get(gca, 'ylim');
    ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
    xlim([0,100])
    set(gca,'XTick',[0:20:100])
    
    subplot(nv,nh,2*nh+2)
    title('\rm Knie')
    hold on
    plot(x,R.Tid(:,iknee),'Color',cols{i},'LineWidth',1)
    axis tight
    yl = get(gca, 'ylim');
    ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
    xlim([0,100])
%     set(gca,'YTick',[-20,0,30,60])
    
    subplot(nv,nh,2*nh+3)
    title('\rm Enkel')
    hold on
    plot(x,R.Tid(:,iankle),'Color',cols{i},'LineWidth',1)
    axis tight
    yl = get(gca, 'ylim');
    ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
    xlim([0,100])
%     set(gca,'YTick',[-90,-50,0,25])
    
    if ~isempty(imtj)
        subplot(nv,nh,2*nh+5)
        title('\rm Midtarsaal gewricht')
        hold on
        plot(x,R.Tid(:,imtj),'Color',cols{i},'LineWidth',1)
        axis tight
        yl = get(gca, 'ylim');
        ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
        xlim([0,100])
        set(gca,'XTick',[0:20:100])
        xlabel('Stapcyclus (%)','Fontsize',label_fontsize);
%         set(gca,'YTick',[-60,-30,0])
    end
    
    subplot(nv,nh,2*nh+6)
    title('\rm Teenbasisgewricht')
    hold on
    plot(x,R.Tid(:,imtp),'Color',cols{i},'LineWidth',1)
    axis tight
    yl = get(gca, 'ylim');
    ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
    xlim([0,100])
    set(gca,'XTick',[0:20:100])
    xlabel('Stapcyclus (%)','Fontsize',label_fontsize);
%     set(gca,'YTick',[-20,-10,0])
    
    %%
%     subplot(nv,nh,2*nh)
%     plot(R.a(:,iSol),'-','Color',cols{i});
%     title('\rm Soleus')
%     ylabel('Activiteit (-)');
%     grid on
% 	set(gca,'XTick',[0:20:100])
%     axis tight
%     yl = get(gca, 'ylim');
%     ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)])
%     xlim([0,100])
        
    %%
    
    subplot(nv,nh,3*nh)
    hold on
    plot(R.GRFs(:,2),'-','Color',cols{i},'LineWidth',1)
    axis tight
    yl = get(gca, 'ylim');
    ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
    xlim([0,100])
    title('\rm Verticale reactiekracht')
    ylabel('F (% gewicht)')
    xlabel('Stapcyclus (%)','Fontsize',label_fontsize);
    
    %%
    
    subplot(nv,nh,4*nh+1)
    hold on
    plot(x,P_joints,'Color',cols{i},'LineWidth',1)
    ylabel('P_{mech} (W/kg)','Fontsize',label_fontsize);
    xlabel('Stapcyclus (%)','Fontsize',label_fontsize);
    title('\rm Voetgewrichten')
    axis tight
    yl = get(gca, 'ylim');
    ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
    xlim([0,100])
    set(gca,'XTick',[0:20:100])
    
    subplot(nv,nh,4*nh+2)
    hold on
    plot(x,P_ankle,'Color',cols{i},'LineWidth',1)
%     ylabel('P_{mech} (W/kg)')
    xlabel('Stapcyclus (%)','Fontsize',label_fontsize);
    title('\rm Enkel')
    axis tight
    yl = get(gca, 'ylim');
    ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
    xlim([0,100])
    set(gca,'XTick',[0:20:100])
    
    if ~isempty(imtj)
        subplot(nv,nh,4*nh+3)
        hold on
        plot(x,P_mtj,'Color',cols{i},'LineWidth',1)
%         ylabel('P_{mech} (W/kg)')
        xlabel('Stapcyclus (%)','Fontsize',label_fontsize);
        title('\rm Midtarsaal gewricht')
        axis tight
        yl = get(gca, 'ylim');
        ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
        xlim([0,100])
        set(gca,'XTick',[0:20:100])
    end
    
    subplot(nv,nh,4*nh+4)
    hold on
    plot(x,P_mtp,'Color',cols{i},'LineWidth',1)
%     ylabel('P_{mech} (W/kg)')
    xlabel('Stapcyclus (%)','Fontsize',label_fontsize);
    title('\rm Teenbasisgewricht')
    axis tight
    yl = get(gca, 'ylim');
    ylim([yl(1)-0.1*abs(sum(yl)),yl(2)+0.1*abs(sum(yl))])
    xlim([0,100])
    set(gca,'XTick',[0:20:100])
%     set(gca,'YTick',[0:20:100])
    
    
end


%%

subplot(nv,nh,4)
hold on
Result_compl_st = load(fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Gefen2002_Gefen2002_Q-30_30_F0_3000_WLv3_ls150.mat']),'R');
R = Result_compl_st.R;
Fs_tib = R.Fs_tib;
l_fa = R.l_fa((R.Qs_mtp(:)==0),:);
plot((l_fa-l_fa(1))*1000,Fs_tib/1000,'--','Color',cols{2},'DisplayName','Intacte voet (soepel)')



Result_compl_noPF_st = load(fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_none_Gefen2002_Q-30_30_F0_3000_WLv3_ls150.mat']),'R');
R = Result_compl_noPF_st.R;
Fs_tib = R.Fs_tib;
l_fa = R.l_fa((R.Qs_mtp(:)==0),:);
plot((l_fa-R.L0)*1000,Fs_tib/1000,':','Color',cols{2},'DisplayName','Zonder plantaire fascia')



Result_stiff_st = load(fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_k300_Q-30_30_F0_3000_WLv3_ls150_sb1_PFx10.mat']),'R');
R = Result_stiff_st.R;
Fs_tib = R.Fs_tib;
l_fa = R.l_fa((R.Qs_mtp(:)==0),:);
plot((l_fa-l_fa(1))*1000,Fs_tib/1000,'--','Color',cols{3},'DisplayName','Intacte voet (stug)')


Result_stiff_noPF_st = load(fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_none_k300_Q-30_30_F0_3000_WLv3_ls150_sb1_PFx10.mat']),'R');
R = Result_stiff_noPF_st.R;
Fs_tib = R.Fs_tib;
l_fa = R.l_fa((R.Qs_mtp(:)==0),:);
plot((l_fa-R.L0)*1000,Fs_tib/1000,':','Color',cols{3},'DisplayName','Zonder plantaire fascia')

title({'Statische simulatie','\rm Drukproef','\rm kadavervoet [2]'})
xlabel('   Horizontale verlenging (mm)','Fontsize',label_fontsize);
ylabel('Verticale kracht (kN)','Fontsize',label_fontsize);
set(gca,'YAxisLocation','right')
xlim([0,9])

%%
lh1=legend('location','northwest');
lhPos = lh1.Position;
lhPos(1) = lhPos(1)-0.01;
lhPos(2) = lhPos(2)-0.173;
set(lh1,'position',lhPos);
% title(lh1,'Statische simulatie')
lh1.Box = 'off';


%%
annotation(gcf,'rectangle',[0.73 0.64 0.22 0.35]);

%%
str = {'\fontsize{8}[1] Data behoort toe aan: A. Falisse, et al.,“Rapid predictive simulations with complex musculoskeletal models (...),”Journal of The Royal Society Interface, vol. 16, no. 157, p. 20190402,2019.',...
    '\fontsize{8}[2] Naar een experiment van: R. F. Ker, et al.,“The spring in the arch of the human foot,”Nature (London), vol. 325, no. 6100,pp. 147–149, 1987.'};
annotation(gcf,'textbox',[0.01,0.01,0.9,0.05],'String',str,'FitBoxToText','on');


%%
h1=gcf;
figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\';
set(h1,'PaperPositionMode','auto')
name = 'one_figure_to_rule_them_all';
print(h1,[figNamePrefix name],'-dpng','-r0')
% print(h1,[figNamePrefix name],'-depsc')

%%
% figure
% pathRefImg = fullfile(pathRepo,'\Figures\ankle_Pothrat.png');
% img_ankle = imread(pathRefImg);
% hold on
% hi1 = image([0,100],flip([-20,25]),img_ankle);
% uistack(hi1,'bottom')
% 
% title('\rm Enkel')
% hold on
% plot(x,R.Qs(:,iankle),'Color',cols{i},'LineWidth',1)
% axis tight











