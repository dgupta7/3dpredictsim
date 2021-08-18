clear
clc
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/VariousFunctions']);
addpath([pathRepo '/PassiveMoments']);
addpath([pathRepo '/FootModel']);

scs = get(0,'ScreenSize');
fsq = [scs(3)/2, scs(4)*0.6];
h1 = figure('Position',[scs(3)/2,140,fsq*0.5]);

resultFiles = {fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_fitted4_Q-20_30_F0_0_WLv3_ls148_mtp1.mat'])
               fullfile([pathRepo '\Results\FootModel\Foot_3D_Fal_s1_mtj_subt1_v5_Natali2010_fitted4_Q-20_30_F0_0_WLv3_ls148_mtp2.mat'])};



load(resultFiles{1},'R');
PlotResults_FootSim(R,[0, 0.4470, 0.7410],h1,1);

load(resultFiles{2},'R');
PlotResults_FootSim(R,[0.4660, 0.6740, 0.1880] 	,h1,1);

xlabel('Passive mtp angle (°)')

lg=legend({'Simulation, R = 9.5 mm','y = 0.987 - 0.000761 x',...
    'Simulation, R = 7.5 mm','y = 0.991 - 0.000513 x'},'Location','southwest');
title(lg,'Metatarsal head radius')


figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\matlab_final\FootSim';
set(h1,'PaperPositionMode','auto')
print(h1,[figNamePrefix '_mtth'],'-dpng','-r0')
print(h1,[figNamePrefix '_mtth'],'-depsc')
    

%% Ker
scs = get(0,'ScreenSize');
set(h,'Position',[scs(3)/2,400,scs(3)/4, scs(3)/4]);
legend({'PF: Gefen2002; mtj: Gefen2002','PF: none; mtj: Gefen2002','PF: none; mtj: Ker1987','PF: Natali2010; mtj: k = 300 Nm/rad'},'Location','northwest','Box','on');
xlim([-0.2,9])
ylim([-0.2,3.4])

figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\matlab_final\FootSim';
set(h,'PaperPositionMode','auto')
print(h,[figNamePrefix '_Ker'],'-dpng','-r0')
print(h,[figNamePrefix '_Ker'],'-depsc')

%% Welte
scs = get(0,'ScreenSize');
set(h,'Position',[scs(3)/2,400,scs(3)/4, scs(3)/4]);

lh=legend({'PF: Gefen2002; mtj: Gefen2002 (-30° mtp dor)','PF: Gefen2002; mtj: Gefen2002 (+30° mtp dor)',...
    'PF: 2*Gefen2002; mtj: fitted6 (-30° mtp dor)','PF: 2*Gefen2002; mtj: fitted6 (+30° mtp dor)',...
    'PF: Natali2010; mtj: k = 300 Nm/rad (-30° mtp dor)','PF: Natali2010; mtj: k = 300 Nm/rad (+30° mtp dor)'},...
    'Location','southoutside');
title(lh,'');
xlim([-0.05,1.1])
ylim([-0.05,1.1])

figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\matlab_final\FootSim';
set(h,'PaperPositionMode','auto')
print(h,[figNamePrefix '_Welte'],'-dpng','-r0')
print(h,[figNamePrefix '_Welte'],'-depsc')


%%

scs = get(0,'ScreenSize');
set(h,'Position',[scs(3)/2,400,scs(3)/4, scs(4)/4]);
lh=legend({'Natali2010','Song2011'},'Location','northwest');

figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\matlab_final\FootSim_PF_GRF';
set(h,'PaperPositionMode','auto')
print(h,[figNamePrefix],'-dpng','-r0')
print(h,[figNamePrefix],'-depsc')
ylabel('F_P_F (N)')


%%

scs = get(0,'ScreenSize');
set(h,'Position',[scs(3)/2,400,scs(3)/4, scs(4)/4]);
lh=legend({'Falisse 2019','Reule 2010','Parr 2012'},'Location','northeast');
title(lh,'Subtalar axis orientation')

figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\matlab_final\FootSim_subt';
set(h,'PaperPositionMode','auto')
print(h,[figNamePrefix],'-dpng','-r0')
print(h,[figNamePrefix],'-depsc')
% ylabel('Angle (°)')

%% Welte
scs = get(0,'ScreenSize');
set(h,'Position',[scs(3)/2,400,scs(3)/4, scs(3)/4]);

lh=legend({'Normal subtalar range (-30° mtp dor)','Normal subtalar range (+30° mtp dor)',...
    '1/2 subtalar range (-30° mtp dor)','1/2 subtalar range (+30° mtp dor)',...
    '1/3 subtalar range (-30° mtp dor)','1/3 subtalar range (+30° mtp dor)'},...
    'Location','southoutside');
title(lh,'Legend');
xlim([-0.05,1.1])
ylim([-0.05,1.1])

figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\matlab_final\FootSim';
set(h,'PaperPositionMode','auto')
print(h,[figNamePrefix '_Welte_subt'],'-dpng','-r0')
print(h,[figNamePrefix '_Welte_subt'],'-depsc')

%%

figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\matlab_final\FootSim';
set(h,'PaperPositionMode','auto')
print(h,[figNamePrefix '_balancing'],'-dpng','-r0')
print(h,[figNamePrefix '_balancing'],'-depsc')

