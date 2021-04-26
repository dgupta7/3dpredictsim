% function [] = PlotResults_3DSim_gaitPhases(ResultsFile)
clear
clc

% ResultsFile = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_k300_MTP_T5_ig24_pp.mat';
ResultsFile = 'D:\school\WTK\thesis\model\3dpredictsim\Results\MidTarsalJoint\Fal_s1_bCst_PF_Natali2010_ls150_MT_nl_signed_lin_MTP_T5_ig24_300_30_pp.mat';

AddCasadiPaths();

pathmain        = pwd;
[pathRepo,~,~]  = fileparts(pathmain);
addpath(genpath(pathRepo));

%% Load external functions
import casadi.*
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB. More details
% about the external function can be found in the documentation.

% Loading external functions.
pathExternalFunctions = [pathRepo,'/ExternalFunctions'];
cd(pathExternalFunctions)

F  = external('F','Leg_3D_Fal_s1_mtj_pp.dll'); 
cd(pathmain);

%%
load(ResultsFile,'R');

for i=1:size(R.GRFs,1)
   GRF(i) = norm(R.GRFs(i,1:3));
end
mg = max(GRF);
grey = [1,1,1]/2;

idx_Qs = [1,2,3,4,5,6,10,11,12,14,16,18,20,22];

th = 3;
i_hs = find(R.GRFs_separate(:,2)>th,1,'first'); % heel strike
i_ho = find(R.GRFs_separate(:,2)>th,1,'last'); % heel off
i_fs = find(R.GRFs_separate(:,5)>th,1,'first'); % forefoot strike
i_fo = find(R.GRFs_separate(:,5)>th,1,'last'); % forefoot off
i_ts = find(R.GRFs_separate(:,8)>th,1,'first'); % toe strike
i_to = find(R.GRFs_separate(:,8)>th,1,'last'); % toe off
ms = round((i_hs+i_to)/2); % mid stance

% idx_ph = [i_hs,i_fs,i_ts,ms,i_ho,i_fo,i_to]; % points of interest
idx_ph = round(linspace(i_hs,i_to,10)); % equally spaced

figure('Position',[100,200,1300,600]);
for i=1:length(idx_ph)
    Qs = R.Qs(idx_ph(i),idx_Qs)'*pi/180;
    Qs(4:6) = Qs(4:6)*180/pi;
    Qdds = zeros(length(Qs),1);
    QsQds = zeros(length(Qs)*2,1);
    QsQds(1:2:end) = Qs;
    res = full(F([QsQds;Qdds]));
    
    leg.x(i,:) = [Qs(4), res([15,18,21])']+(i-1)*0.25;
    leg.y(i,:) = [Qs(5), res([16,19,22])'];
    
    hindfoot.x(i,:) = [res([21,27])', res([36,33])'+res(24),res(21)]+(i-1)*0.25;
    hindfoot.y(i,:) = [res([22,28])', res([37,34])'+res(25),res(22)];
    
    forefoot.x(i,:) = [res([27,30])', res([42,39])'+res(27),res(27)]+(i-1)*0.25;
    forefoot.y(i,:) = [res([28,31])', res([43,40])'+res(28),res(28)];
    
    toes.x(i,:) = [res(30)', res([48,45])'+res(30),res(30)]+(i-1)*0.25;
    toes.y(i,:) = [res(31)', res([49,46])'+res(31),res(31)];

    Qs(end) = 0; % mtp = 0
    QsQds(1:2:end) = Qs;
    res2 = full(F([QsQds;Qdds]));
    toes.x0(i,:) = [res2(30)', res2([48,45])'+res2(30),res2(30)]+(i-1)*0.25;
    toes.y0(i,:) = [res2(31)', res2([49,46])'+res2(31),res2(31)];
    
    Qs(end-1) = 0; % mtj = 0
    QsQds(1:2:end) = Qs;
    res3 = full(F([QsQds;Qdds]));
    forefoot.x0(i,:) = [res3([27,30])', res3([42,39])'+res3(27),res3(27)]+(i-1)*0.25;
    forefoot.y0(i,:) = [res3([28,31])', res3([43,40])'+res3(28),res3(28)];
    
    Qs(end-2) = 0; % subt = 0
    Qs(end-3) = 0; % ankle = 0
    QsQds(1:2:end) = Qs;
    res4 = full(F([QsQds;Qdds]));
    hindfoot.x0(i,:) = [res4([21,27])', res4([36,33])'+res4(24),res4(21)]+(i-1)*0.25;
    hindfoot.y0(i,:) = [res4([22,28])', res4([37,34])'+res4(25),res4(22)];
    
    subplot(4,1,[1:3])
    hold on
    plot(leg.x(i,:),leg.y(i,:),'k')
    plot(leg.x(i,:),leg.y(i,:),'.k')
    plot(hindfoot.x(i,:),hindfoot.y(i,:),'k')
    plot(hindfoot.x(i,end),hindfoot.y(i,end),'.k')
    plot(forefoot.x(i,:),forefoot.y(i,:),'k')
    plot(forefoot.x(i,end),forefoot.y(i,end),'.k')
    plot(toes.x(i,:),toes.y(i,:),'k')
    axis equal
    
    leg_p.x(i,:) = [(leg.x(i,end-1)+4*leg.x(i,end))/5, leg.x(i,end)];
    leg_p.y(i,:) = [(leg.y(i,end-1)+4*leg.y(i,end))/5, leg.y(i,end)];
    
    grfs(i) = norm(R.GRFs(idx_ph(i),1:3));
    
    grf.x(i,:) = [R.COPR(idx_ph(i),1),...
        R.COPR(idx_ph(i),1) - R.GRFs(idx_ph(i),1)/grfs(i)/10] +(i-1)*0.25;
    grf.y(i,:) = [R.COPR(idx_ph(i),2),...
        R.COPR(idx_ph(i),2) - R.GRFs(idx_ph(i),2)/grfs(i)/10];
    
    col(1) = interp1([0,0.25,0.5,0.75,1],[0,0,1,1,1],grfs(i)/mg);
    col(2) = interp1([0,0.25,0.5,0.75,1],[0,1,1,0.5,0],grfs(i)/mg);
    col(3) = interp1([0,0.25,0.5,0.75,1],[1,0,0,0,0],grfs(i)/mg);
    
    cols(i,:) = col';
    
    
    subplot(4,10,31:39)
    hold on
    plot(hindfoot.x0(i,:),hindfoot.y0(i,:),'Color',grey)
    plot(forefoot.x0(i,:),forefoot.y0(i,:),'Color',grey)
    plot(toes.x0(i,:),toes.y0(i,:),'Color',grey)
    
    plot(leg_p.x(i,:),leg_p.y(i,:),'k')
    plot(leg_p.x(i,end),leg_p.y(i,end),'.k')
    plot(hindfoot.x(i,:),hindfoot.y(i,:),'k')
    plot(hindfoot.x(i,end),hindfoot.y(i,end),'.k')
    plot(forefoot.x(i,:),forefoot.y(i,:),'k')
    plot(forefoot.x(i,end),forefoot.y(i,end),'.k')
    plot(toes.x(i,:),toes.y(i,:),'k')
    plot(grf.x(i,:),grf.y(i,:),'Color',col)
    plot(grf.x(i,1),grf.y(i,1),'^','Color',col)
    axis equal
    
end

subplot(4,1,[1:3])
line(get(gca, 'xlim'),[0,0],'color','k','LineStyle','--')

subplot(4,10,31:39)
line(get(gca, 'xlim'),[0,0],'color','k','LineStyle','--')

subplot(4,10,40)
hold on
yaxis right
for i=1:length(idx_ph)
    plot(0,grfs(i),'d','Color',cols(i,:),'MarkerFace',cols(i,:))
end
ylabel('GRF (%BW)')




% end