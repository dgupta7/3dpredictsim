function [] = PlotResults_3DSim(ResultsFile,Cs,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

load(ResultsFile);
% create figure handles
if ~isempty(varargin)
    h = varargin{1};
    figure(h);
    tab1 = h.Parent.Children(1).Children(1).Children(1);
    tab2 = h.Parent.Children(1).Children(1).Children(2);
    tab3 = h.Parent.Children(1).Children(1).Children(3);
    tab4 = h.Parent.Children(1).Children(1).Children(4);
    tab5 = h.Parent.Children(1).Children(1).Children(5);
else
    h = figure();
    hTabGroup = uitabgroup;
    tab1 = uitab(hTabGroup, 'Title', 'MainInfo');
    tab2 = uitab(hTabGroup, 'Title', 'ExoT');
    tab3 = uitab(hTabGroup, 'Title', 'CalfM');
    tab4 = uitab(hTabGroup, 'Title', 'Kinematics');
    tab5 = uitab(hTabGroup, 'Title', 'Kinetics');
end


%% Plot general information
axes('parent', tab1);

% Plot COT as a function of exo assitance
subplot(2,2,1); hold on;
plot(R.Sopt.ExoScale,R.COT,'o','Color',Cs,'MarkerFaceColor',Cs);
ylabel('COT');  xlabel('Level assistance');
title('Cost of Transport');

% Plot stride frequency
subplot(2,2,2); hold on;
dt = R.t(end);
plot(R.Sopt.ExoScale,1./dt,'o','Color',Cs,'MarkerFaceColor',Cs); 
ylabel('Stride Frequency'); xlabel('Level assistance');
title('Stride Frequency');

% 
subplot(2,2,3);  hold on;
plot(R.T_exo(:,1),'-','Color',Cs);
ylabel('Exo Moment [Nm]');  xlabel('% stride');
title('Left');
% text(20,-10,['Time shift: ' num2str(round(R.dt_exoShift,3)) ' s']);


subplot(2,2,4); hold on;
plot(R.T_exo(:,2),'-','Color',Cs); 
ylabel('Exo Moment [Nm]');  xlabel('% stride');
title('Right');


%% Plot Torque information
axes('parent', tab2);

subplot(3,2,1)
plot(R.T_exo(:,1),'-','Color',Cs); hold on;
ylabel('Exo Moment [Nm]');  xlabel('% stride');
title('Left');

subplot(3,2,2);
plot(R.T_exo(:,2),'-','Color',Cs); hold on;
ylabel('Exo Moment [Nm]'); xlabel('% stride');
title('Right');

subplot(3,2,3)
plot(R.Tid(:,strcmp(R.colheaders.joints,'ankle_angle_l')),'-','Color',Cs); hold on;
ylabel('Ankle moment [Nm]'); xlabel('% stride');
title('Left');

subplot(3,2,4)
plot(R.Tid(:,strcmp(R.colheaders.joints,'ankle_angle_r')),'-','Color',Cs); hold on;
ylabel('Ankle moment [Nm]'); xlabel('% stride');
title('Right');

subplot(3,2,5)
plot(R.Tid(:,strcmp(R.colheaders.joints,'ankle_angle_l'))-R.T_exo(:,1),'-','Color',Cs); hold on;
ylabel('Biological ankle moment [Nm]'); xlabel('% stride');
title('Left');

subplot(3,2,6)
plot(R.Tid(:,strcmp(R.colheaders.joints,'ankle_angle_r'))-R.T_exo(:,2),'-','Color',Cs); hold on;
ylabel('Biological ankle moment [Nm]'); xlabel('% stride');
title('Left');

ax =[];
ax2 = [];
for i=1:3
   ax(i) = subplot(3,2,i*2-1); 
   ax2(i) = subplot(3,2,i*2); 
end
linkaxes(ax,'x');
linkaxes(ax2,'x');

%% Plot ankle muscle energetics
axes('parent', tab3);
iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
iGas = find(strcmp(R.colheaders.muscles,'lat_gas_r'));

subplot(5,2,1)
plot(R.a(:,iSol),'-','Color',Cs); hold on; title('Soleus');
xlabel('% stride'); ylabel('activity');

subplot(5,2,2)
plot(R.a(:,iGas),'-','Color',Cs); hold on; title('Gastrocnemius');
xlabel('% stride'); ylabel('activity');

subplot(5,2,3)
plot(R.MetabB.Etot(:,iSol),'-','Color',Cs); hold on;
xlabel('% stride'); ylabel('Muscle metab power');

subplot(5,2,4)
plot(R.MetabB.Etot(:,iGas),'-','Color',Cs); hold on;
xlabel('% stride'); ylabel('Muscle metab power (W)');

subplot(5,2,5)
plot(R.lMtilde(:,iSol),'-','Color',Cs); hold on;
xlabel('% stride'); ylabel('Norm fiber length');

subplot(5,2,6)
plot(R.lMtilde(:,iGas),'-','Color',Cs); hold on;
xlabel('% stride'); ylabel('Norm fiber length');

subplot(5,2,7)
plot(R.MetabB.Wdot(:,iSol),'-','Color',Cs); hold on;
xlabel('% stride'); ylabel('Wdot');

subplot(5,2,8)
plot(R.MetabB.Wdot(:,iGas),'-','Color',Cs); hold on;
xlabel('% stride'); ylabel('Wdot');

subplot(5,2,9)
plot(R.FT(:,iSol),'-','Color',Cs); hold on;
xlabel('% stride'); ylabel('Norm muscle force');

subplot(5,2,10)
plot(R.FT(:,iGas),'-','Color',Cs); hold on;
xlabel('% stride'); ylabel('Norm muscle force');

ax =[];
ax2 = [];
for i=1:5
   ax(i) = subplot(5,2,i*2-1); 
   ax2(i) = subplot(5,2,i*2); 
end
linkaxes(ax,'x');
linkaxes(ax2,'x');



end

