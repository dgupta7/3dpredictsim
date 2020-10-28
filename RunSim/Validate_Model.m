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

[ccc1] = ModelValidation(ResultsFile1, pathData);
[ccc2] = ModelValidation(ResultsFile2, pathData);

%% Plot

h = figure();
h.Name = 'ModelValidation';
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
    stem(ccc1.kinematics.shift(:,i),ccc1.kinematics.ccc(:,i),'linewidth',line_linewidth,'DisplayName','b cst')
    hold on
    stem(ccc2.kinematics.shift(:,i),ccc2.kinematics.ccc(:,i),'linewidth',line_linewidth,'DisplayName','alpha cst')
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
    stem(ccc1.kinetics.shift(:,i),ccc1.kinetics.ccc(:,i),'linewidth',line_linewidth,'DisplayName','b cst')
    hold on
    stem(ccc2.kinetics.shift(:,i),ccc2.kinetics.ccc(:,i),'linewidth',line_linewidth,'DisplayName','alpha cst')
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











saveas(h,fullfile(pathFig,'FigureModelValidation.fig'));


%%
addpath([pathRepo,'/Plots']);
hh = figure(); 	% new figure with handle
PlotResults_3DSim(ResultsFile1,[1 0 0],'b cst',hh);
PlotResults_3DSim(ResultsFile2,[0 0 1],'alpha cst',hh);





