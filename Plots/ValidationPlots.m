function [] = ValidationPlots(pathData,ResultsFile1,varargin)

nr = 1 + length(varargin);

names={nr};
ccc={nr};
NLSE={nr};

[~,names{1},~] = fileparts(ResultsFile1);
ccc{1} = ModelValidation_CrossCorrelationCoefficient(ResultsFile1, pathData);
NLSE{1} = ModelValidation_NLSError(ResultsFile1, pathData);

for i=2:nr
    [~,names{i},~] = fileparts(varargin{i-1});
    ccc{i} = ModelValidation_CrossCorrelationCoefficient(varargin{i-1}, pathData);
    NLSE{i} = ModelValidation_NLSError(varargin{i-1}, pathData);
end



h = figure();
set(h,'Position',[82         151        1497         827]);
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
for i = 1:length(ccc{1}.kinematics.max)
    subplot(3,6,i)
    plot(ccc{1}.kinematics.shift(:,i),ccc{1}.kinematics.ccc(:,i),'linewidth',line_linewidth,'DisplayName',names{1})
    hold on
    for j=2:nr
        plot(ccc{j}.kinematics.shift(:,i),ccc{j}.kinematics.ccc(:,i),'linewidth',line_linewidth,'DisplayName',names{j})
    end
    grid on
    title(ccc{1}.joints{i},'Fontsize',label_fontsize)
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
% lhPos(1) = lhPos(1)+0.2;
set(lh,'position',lhPos);


% plot kinetics
axes('parent', tab2);  
for i = 1:size(ccc{1}.kinetics.ccc,2)
    subplot(3,6,i)
    plot(ccc{1}.kinetics.shift(:,i),ccc{1}.kinetics.ccc(:,i),'linewidth',line_linewidth,'DisplayName',names{1})
    hold on
    for j=2:nr
        plot(ccc{j}.kinetics.shift(:,i),ccc{j}.kinetics.ccc(:,i),'linewidth',line_linewidth,'DisplayName',names{j})
    end
    grid on
    title(ccc{1}.joints{i},'Fontsize',label_fontsize)
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
% lhPos(1) = lhPos(1)+0.2;
set(lh,'position',lhPos);

% saveas(h,fullfile(pathFig,'FigureModelValidation.fig'));


%% Nonlinear least squares error
% of the simulation results and the mean measurement data, weighted with
% the inverse of the variance



h2 = figure();
set(h2,'Position',[82         151        1497         827]);
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

    
% plot kinematics
axes('parent', h2tab1);  
for i = 1:length(ccc{1}.kinematics.max)
    subplot(3,6,i)
%     plot(0,NLSE1.kinematics(i),'*','DisplayName','b cst')
    plot(0,1,'*','DisplayName',names{1})
    hold on
     for j=2:nr
%     plot(0,NLSE2.kinematics(i),'*','DisplayName','alpha cst')
        plot(0,NLSE{j}.kinematics(i)/NLSE{1}.kinematics(i),'*','DisplayName',names{j})
     end
    grid on
    title(NLSE{1}.joints{i},'Fontsize',label_fontsize)
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
for i = 1:length(ccc{1}.kinetics.max)
    subplot(3,6,i)
%     plot(0,NLSE1.kinetics(i),'*','DisplayName','b cst')
    plot(0,1,'*','DisplayName',names{1})
    hold on
    for j=2:nr
%     plot(0,NLSE2.kinetics(i),'*','DisplayName','alpha cst')
        plot(0,NLSE{j}.kinetics(i)/NLSE{1}.kinetics(i),'*','DisplayName',names{j})
    end
    grid on
    title(NLSE{1}.joints{i},'Fontsize',label_fontsize)
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
% lhPos(1) = lhPos(1)+0.2;
set(lh,'position',lhPos);


%  if boolFirst && md
%             iSol_data = find(strcmp(Dat.Normal.EMGheaders,'soleus_r'));
%             iGas_data = find(strcmp(Dat.Normal.EMGheaders,'gas_lat_r'));
%             sol_act_normal = [Dat.Normal.gc.lowEMG_mean(ceil(end/2):end,iSol_data); Dat.Normal.gc.lowEMG_mean(1:ceil(end/2)-1,iSol_data)];
%             sol_act_active = [Dat.Active.gc.lowEMG_mean(ceil(end/2):end,iSol_data); Dat.Active.gc.lowEMG_mean(1:ceil(end/2)-1,iSol_data)];
%             gas_act_normal = [Dat.Normal.gc.lowEMG_mean(ceil(end/2):end,iGas_data); Dat.Normal.gc.lowEMG_mean(1:ceil(end/2)-1,iGas_data)];
%             gas_act_active = [Dat.Active.gc.lowEMG_mean(ceil(end/2):end,iGas_data); Dat.Active.gc.lowEMG_mean(1:ceil(end/2)-1,iGas_data)];
%                         
%     end

%     if boolFirst && md
%         plot(sol_act_normal,'-','Color',darkgray,'DisplayName','Data normal shoes');
%         plot(sol_act_active,'-','Color',lightgray,'DisplayName','Data exo active');
%     end

end