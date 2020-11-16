clear 
close all
clc


%% Run

[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/VariousFunctions']);
pathResults = [pathRepo,'/Results'];
pathData = [pathRepo,'/ExperimentalData','/ExperimentalData.mat'];
pathFig = [pathRepo,'/Figures/Test_Lars'];


ResultsFile1 = fullfile(pathResults,'Test_Lars','Test1_bCst_no_pp.mat');
ResultsFile2 = fullfile(pathResults,'Test_Lars','Test1_alphaCst_no_pp.mat');

ResultsFile3 = fullfile(pathResults,'Test_Lars','Test1_bCst_pas_pp.mat');
ResultsFile4 = fullfile(pathResults,'Test_Lars','Test1_alphaCst_pas_pp.mat');

ResultsFile5 = fullfile(pathResults,'Test_Lars','Test1_bCst_exo_pp.mat');
ResultsFile6 = fullfile(pathResults,'Test_Lars','Test1_alphaCst_exo_pp.mat');

%% Cross-correlation coefficient
% express similarity in shape of the simulation results and the mean measurement data

[ccc1] = ModelValidation_CrossCorrelationCoefficient(ResultsFile1, pathData);
[ccc2] = ModelValidation_CrossCorrelationCoefficient(ResultsFile2, pathData);

h = figure();
h.Name = 'ModelValidation: Cross-correlation coefficient';
hTabGroup = uitabgroup;
tab1 = uitab(hTabGroup, 'Title', 'Kinematics');
tab2 = uitab(hTabGroup, 'Title', 'Kinetics');
% tab3 = uitab(hTabGroup, 'Title', 'COT');
% tab4 = uitab(hTabGroup, 'Title', 'ExoInfo');
% tab5 = uitab(hTabGroup, 'Title', 'CalfM');
% tab6 = uitab(hTabGroup, 'Title', 'Ground reaction force');
% tab7 = uitab(hTabGroup, 'Title', 'Objective Function');
% tab8 = uitab(hTabGroup, 'Title', 'Ankle detailed');
% tab9 = uitab(hTabGroup, 'Title', 'SpatioTemporal');
set(h,'Color','w');
label_fontsize  = 12;
line_linewidth  = 0.5;
    
% plot kinematics
axes('parent', tab1);  
for i = 1:length(ccc1.kinematics.max)
    subplot(3,6,i)
    plot(ccc1.kinematics.shift(:,i),ccc1.kinematics.ccc(:,i),'linewidth',line_linewidth,'DisplayName','b cst')
    hold on
    plot(ccc2.kinematics.shift(:,i),ccc2.kinematics.ccc(:,i),'linewidth',line_linewidth,'DisplayName','alpha cst')
    grid on
    title(ccc1.joints{i},'Fontsize',label_fontsize)
    if i == 1 || i == 7 ||i == 13
        ylabel('cross-correlation coefficient','Fontsize',label_fontsize);
    end
    if i > 12     
        xlabel('shift (% gait cycle)','Fontsize',label_fontsize);
    end
end
lh=legend('-DynamicLegend','location','east');
lh.Interpreter = 'none';
lhPos = lh.Position;
lhPos(1) = lhPos(1)+0.2;
set(lh,'position',lhPos);


% plot kinetics
axes('parent', tab2);  
for i = 1:length(ccc1.kinetics.max)
    subplot(3,6,i)
    plot(ccc1.kinetics.shift(:,i),ccc1.kinetics.ccc(:,i),'linewidth',line_linewidth,'DisplayName','b cst')
    hold on
    plot(ccc2.kinetics.shift(:,i),ccc2.kinetics.ccc(:,i),'linewidth',line_linewidth,'DisplayName','alpha cst')
    grid on
    title(ccc1.joints{i},'Fontsize',label_fontsize)
    if i == 1 || i == 7 ||i == 13
        ylabel('cross-correlation coefficient','Fontsize',label_fontsize);
    end
    if i > 12     
        xlabel('shift (% gait cycle)','Fontsize',label_fontsize);
    end
end
lh=legend('-DynamicLegend','location','east');
lh.Interpreter = 'none';
lhPos = lh.Position;
lhPos(1) = lhPos(1)+0.2;
set(lh,'position',lhPos);

% saveas(h,fullfile(pathFig,'FigureModelValidation.fig'));


%% Nonlinear least squares error
% of the simulation results and the mean measurement data, weighted with
% the inverse of the variance

[NLSE1] = ModelValidation_NLSError(ResultsFile1, pathData);
[NLSE2] = ModelValidation_NLSError(ResultsFile2, pathData);

h2 = figure();
h2.Name = 'ModelValidation: Nonlinear least squares error';
hTabGroup = uitabgroup;
h2tab1 = uitab(hTabGroup, 'Title', 'Kinematics');
h2tab2 = uitab(hTabGroup, 'Title', 'Kinetics');
% tab3 = uitab(hTabGroup, 'Title', 'COT');
% tab4 = uitab(hTabGroup, 'Title', 'ExoInfo');
% tab5 = uitab(hTabGroup, 'Title', 'CalfM');
% tab6 = uitab(hTabGroup, 'Title', 'Ground reaction force');
% tab7 = uitab(hTabGroup, 'Title', 'Objective Function');
% tab8 = uitab(hTabGroup, 'Title', 'Ankle detailed');
% tab9 = uitab(hTabGroup, 'Title', 'SpatioTemporal');
set(h2,'Color','w');
label_fontsize  = 12;
line_linewidth  = 0.5;
    
% plot kinematics
axes('parent', h2tab1);  
for i = 1:length(ccc1.kinematics.max)
    subplot(3,6,i)
%     plot(0,NLSE1.kinematics(i),'*','DisplayName','b cst')
    plot(0,1,'*','DisplayName','b cst')
    hold on
%     plot(0,NLSE2.kinematics(i),'*','DisplayName','alpha cst')
    plot(0,NLSE2.kinematics(i)/NLSE1.kinematics(i),'*','DisplayName','alpha cst')
    grid on
    title(NLSE1.joints{i},'Fontsize',label_fontsize)
    if i == 1 || i == 7 ||i == 13
        ylabel('NLSE','Fontsize',label_fontsize);
    end
%     if i > 12     
%         xlabel('shift (% gait cycle)','Fontsize',label_fontsize);
%     end
end
lh2=legend('-DynamicLegend','location','east');
lh2.Interpreter = 'none';
lhPos = lh2.Position;
lhPos(1) = lhPos(1)+0.2;
set(lh2,'position',lhPos);


% plot kinetics
axes('parent', h2tab2);  
for i = 1:length(ccc1.kinetics.max)
    subplot(3,6,i)
%     plot(0,NLSE1.kinetics(i),'*','DisplayName','b cst')
    plot(0,1,'*','DisplayName','b cst')
    hold on
%     plot(0,NLSE2.kinetics(i),'*','DisplayName','alpha cst')
    plot(0,NLSE2.kinetics(i)/NLSE1.kinetics(i),'*','DisplayName','alpha cst')
    grid on
    title(NLSE1.joints{i},'Fontsize',label_fontsize)
    if i == 1 || i == 7 ||i == 13
        ylabel('NLSE','Fontsize',label_fontsize);
    end
%     if i > 12     
%         xlabel('shift (% gait cycle)','Fontsize',label_fontsize);
%     end
end
lh=legend('-DynamicLegend','location','east');
lh.Interpreter = 'none';
lhPos = lh.Position;
lhPos(1) = lhPos(1)+0.2;
set(lh,'position',lhPos);

% saveas(h2,fullfile(pathFig,'FigureModelValidation.fig'));






%% Make general plots for both results
addpath([pathRepo,'/Plots']);
hh1 = figure(); 	% new figure with handle
PlotResults_3DSim(ResultsFile1,[1 0 0],'b cst (no exo)',hh1);
PlotResults_3DSim(ResultsFile2,[0 0 1],'alpha cst (no exo)',hh1);

hh2 = figure(); 	% new figure with handle
PlotResults_3DSim(ResultsFile3,[1 0 0],'b cst (pas)',hh2);
PlotResults_3DSim(ResultsFile4,[0 0 1],'alpha cst (pas)',hh2);

hh3 = figure(); 	% new figure with handle
PlotResults_3DSim(ResultsFile5,[1 0 0],'b cst (exo)',hh3);
PlotResults_3DSim(ResultsFile6,[0 0 1],'alpha cst (exo)',hh3);


%% Plot absolute and relative difference of simulation wrt mean measurement

% PlotResultsComparison_3DSim(ResultsFile1,ResultsFile2)
% PlotResultsComparison_3DSim(ResultsFile3,ResultsFile4)
% PlotResultsComparison_3DSim(ResultsFile5,ResultsFile6)





