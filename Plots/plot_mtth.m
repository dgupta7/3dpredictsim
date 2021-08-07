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
set(h,'Position',[scs(3)/2,400,scs(3)/4, scs(3)/4]);
legend({'PF: Gefen2002; mtj: Gefen2002','PF: none; mtj: Gefen2002','PF: none; mtj: Ker1987'},'Location','northwest','Box','off');
xlim([-0.2,9])
ylim([-0.2,3.4])

figNamePrefix = 'D:\OneDrive\WTK\thesis\figuren\matlab_final\FootSim';
set(h,'PaperPositionMode','auto')
print(h,[figNamePrefix '_Ker'],'-dpng','-r0')
print(h,[figNamePrefix '_Ker'],'-depsc')
