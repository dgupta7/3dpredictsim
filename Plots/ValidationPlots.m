function [] = ValidationPlots(ResultsFile,Cs,LegName,varargin)

set(0,'defaultTextInterpreter','none');

if ~exist('LegName','var')
    [~,filename,~]= fileparts(ResultsFile);
    LegName = filename;
end

if exist(ResultsFile,'file')
    load(ResultsFile,'R');
    
    
    if ~isfield(R,'CrossCorrelation') || isempty(R.CrossCorrelation)
        if strcmp(R.S.subject,'Fal_s1') || strcmp(R.S.subject,'subject1')
            pathData = 'D:\school\WTK\thesis\model\3dpredictsim\Data\Fal_s1.mat';
        else
            % load data Pog_s1 from struct saved during ...\Analyze_ExoData\Batch\BatchScript_LatexReport.m
            pathData = 'D:\school\WTK\thesis\model\3dpredictsim\Data\Pog_s1.mat';
        end
        R.CrossCorrelation = ModelValidation_CrossCorrelationCoefficient(ResultsFile, pathData);
        save(ResultsFile,'R');
    end

    boolFirst = 1;
    
    % create figure handles
    if ~isempty(varargin) && ~isempty(varargin{1})
        if ~isempty(varargin{1}.Name)
            h = varargin{1};
            tab1 = h.Parent.Children(1).Children(1).Children(1);
            tab2 = h.Parent.Children(1).Children(1).Children(2);
            tab3 = h.Parent.Children(1).Children(1).Children(3);
            tab4 = h.Parent.Children(1).Children(1).Children(4);
            tab5 = h.Parent.Children(1).Children(1).Children(5);
            boolFirst = 0;
        else
            h = varargin{1};
            hTabGroup = uitabgroup;
            tab1 = uitab(hTabGroup, 'Title', 'Kinematics');
            tab2 = uitab(hTabGroup, 'Title', 'Kinetics');
            tab3 = uitab(hTabGroup, 'Title', 'Muscle activity');
            tab4 = uitab(hTabGroup, 'Title', 'Ankle, Soleus');
            tab5 = uitab(hTabGroup, 'Title', 'Knee, Hip');
            h.Name = 'ModelValidation: Cross-correlation coefficient';
            set(h,'Color','w');
        end
    else
        h = figure();
        h.Name = 'ModelValidation: Cross-correlation coefficient';
        hTabGroup = uitabgroup;
        tab1 = uitab(hTabGroup, 'Title', 'Kinematics');
        tab2 = uitab(hTabGroup, 'Title', 'Kinetics');
        tab3 = uitab(hTabGroup, 'Title', 'Muscle activity');
        tab4 = uitab(hTabGroup, 'Title', 'Ankle, Soleus');
        tab5 = uitab(hTabGroup, 'Title', 'Knee, Hip');
        set(h,'Color','w');
    end



    label_fontsize  = 12;
    line_linewidth  = 0.5;
    
    %% plot kinematics
    axes('parent', tab1);

    for i = 1:length(R.CrossCorrelation.kinematics.max)
        subplot(3,6,i)
        hold on
        grid on
        plot(R.CrossCorrelation.kinematics.shift(:,i),R.CrossCorrelation.kinematics.ccc(:,i),'color',Cs,'linewidth',line_linewidth,'DisplayName',LegName)
        if boolFirst
            title(R.CrossCorrelation.joints{i},'Fontsize',label_fontsize)
            if i == 1 || i == 7 ||i == 13
                ylabel('cross-correlation coefficient','Fontsize',label_fontsize);
            end
            if i > 12     
                xlabel('shift (% gait cycle)','Fontsize',label_fontsize);
            end
        end
    end
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
        lhPos = lh.Position;
        % lhPos(1) = lhPos(1)+0.2;
        set(lh,'position',lhPos);
    end

    %% plot kinetics
    axes('parent', tab2);  
    for i = 1:size(R.CrossCorrelation.kinetics.ccc,2)
        subplot(3,6,i)
        hold on
        grid on
        plot(R.CrossCorrelation.kinetics.shift(:,i),R.CrossCorrelation.kinetics.ccc(:,i),'color',Cs,'linewidth',line_linewidth,'DisplayName',LegName)
        if boolFirst
            title(R.CrossCorrelation.joints{i},'Fontsize',label_fontsize)
            if i == 1 || i == 7 ||i == 13
                ylabel('cross-correlation coefficient','Fontsize',label_fontsize);
            end
            if i > 12     
                xlabel('shift (% gait cycle)','Fontsize',label_fontsize);
            end
        end
    end
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
        lhPos = lh.Position;
        % lhPos(1) = lhPos(1)+0.2;
        set(lh,'position',lhPos);
    end

    %% plot muscles
    axes('parent', tab3);
    for i = 1:size(R.CrossCorrelation.muscles.ccc,2)
        subplot(3,3,i)
        hold on
        grid on
        plot(R.CrossCorrelation.muscles.shift(:,i),R.CrossCorrelation.muscles.ccc(:,i),'color',Cs,'linewidth',line_linewidth,'DisplayName',LegName)
        if boolFirst
            title(R.CrossCorrelation.muscles.names{i},'Fontsize',label_fontsize)
            if i == 1 || i == 4 || i == 7
                ylabel('cross-correlation coefficient','Fontsize',label_fontsize);
            end
            if i == 7 || i == 8 || i == 6
                xlabel('shift (% gait cycle)','Fontsize',label_fontsize);
            end
        end
    end
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
        lhPos = lh.Position;
        lhPos(1) = lhPos(1)+0.2;
        set(lh,'position',lhPos);
    end

%     %% interpretation
%     
%     iankle = find(strcmp(R.CrossCorrelation.joints,'Ankle R'));
%     iknee = find(strcmp(R.CrossCorrelation.joints,'Knee R'));
%     ihipa = find(strcmp(R.CrossCorrelation.joints,'Hip adduction R'));
%     ihipf = find(strcmp(R.CrossCorrelation.joints,'Hip flexion R'));
%     isol = find(strcmp(R.CrossCorrelation.muscles.names,'Soleus'));
%     
%     
%     has_tmt_unlocked =  isfield(R.S,'tmt') && ~isempty(R.S.tmt) && R.S.tmt && isfield(R.S,'tmt_locked') && ~isempty(R.S.tmt_locked) && ~R.S.tmt_locked;
%     has_WL = isfield(R.S,'Windlass') && ~isempty(R.S.Windlass) && R.S.Windlass ~= 0;
%     has_mtj = isfield(R.S,'mtj') && ~isempty(R.S.mtj) && R.S.mtj;
%     
%     % ankle and soleus
%     axes('parent', tab4);
% 
%     if has_tmt_unlocked || has_mtj
%         if has_WL
%             mrk = 'o';
%         else
%             mrk = '*';
%             R.S.cWL = 0;
%         end
%         xpar1 = [R.S.kTMT,R.S.cWL,R.S.MuscModelAsmp,R.S.kTMT,R.S.cWL];
%         xpar2 = [R.S.kTMT,R.S.cWL,R.S.MuscModelAsmp,R.S.MuscModelAsmp,R.S.cWL];
%         xpar = [xpar1,xpar1,xpar2,xpar2];
%     end
% 
%     ypar1 = [[1,1,1]*R.CrossCorrelation.kinematics.max(1,iankle), [1,1]*R.CrossCorrelation.muscles.max(1,isol)];
%     ypar2 = [[1,1,1]*R.CrossCorrelation.kinematics.shiftAtMax(1,iankle), [1,1]*R.CrossCorrelation.muscles.shiftAtMax(1,isol)];
%     ypar3 = [[1,1,1]*R.CrossCorrelation.kinetics.max(1,iankle), [1,1]*R.CrossCorrelation.muscles.max(1,isol)];
%     ypar4 = [[1,1,1]*R.CrossCorrelation.kinetics.shiftAtMax(1,iankle), [1,1]*R.CrossCorrelation.muscles.shiftAtMax(1,isol)];
% 
%     ypar = [ypar1,ypar2,ypar3,ypar4];
% 
%     for i=1:19
%         subplot(4,5,i)
%         hold on
%         grid on
%         if has_tmt_unlocked || has_mtj
%             plot(xpar(i),ypar(i),mrk,'Color',Cs,'DisplayName',LegName)
%         else
%             line(get(gca, 'xlim'),[1,1]*ypar(i),'color',Cs,'DisplayName',LegName)
%         end
% 
%         if boolFirst
% 
%             if i > 0 && i < 4
%                 title('Ankle angle','Fontsize',label_fontsize)
%             end
%             if i > 10 && i < 14
%                 title('Ankle torque','Fontsize',label_fontsize)
%             end
%             if (i > 3 && i < 6) || i == 14
%                 title('Soleus activity','Fontsize',label_fontsize)
%             end
%             if i == 1 || i == 11
%                 ylabel('Max R','Fontsize',label_fontsize);
%             end
%             if i == 6 || i == 16
%                 ylabel('R shift (%)','Fontsize',label_fontsize);
%             end
%             if i == 6 || i == 9 || i == 16
%                 xlabel('k tmt (Nm/rad)','Fontsize',label_fontsize);
%             end
%             if i == 7 || i == 10 || i == 17
%                 xlabel('cWL (/°)','Fontsize',label_fontsize);
%             end
%             if i == 8 || i == 19 || i == 18
%                 xlabel('h cst      alpha cst','Fontsize',label_fontsize); 
%             end
%             if i == 15
%                 lh=legend('-DynamicLegend','location','northwest');
%                 lh.Interpreter = 'none';
% 
%             end
%         end
% 
%     end
%     
% 
%     % knee and hip
%     axes('parent', tab5);
%     
%     if has_tmt_unlocked || has_mtj
%         if has_WL
%             mrk = 'o';
%         else
%             mrk = '*';
%             R.S.cWL = 0;
%         end
%         xpar = [R.S.kTMT,R.S.cWL];
%     end
% 
%     ypar = [R.CrossCorrelation.kinematics.max(1,:);
%             R.CrossCorrelation.kinematics.shiftAtMax(1,:);
%             R.CrossCorrelation.kinetics.max(1,:);
%             R.CrossCorrelation.kinetics.shiftAtMax(1,:)];
%     iis = [iknee, ihipf, ihipa];
%     
%     
%     for i=1:24
%         tmp1 = mod(i-1,2)+1;
%         tmp2 = ceil(i/6);
%         tmp3 = ceil((mod(i-1,6)+1)/2);
%         
%         subplot(4,6,i)
%         hold on
%         grid on
%         if has_tmt_unlocked || has_mtj
%             plot(xpar(tmp1),ypar(tmp2,iis(tmp3)),mrk,'Color',Cs,'DisplayName',LegName)
%         else
%             line(get(gca, 'xlim'),[1,1]*ypar(tmp2,iis(tmp3)),'color',Cs,'DisplayName',LegName)
%         end
% 
%         if boolFirst
% 
%             if i > 0 && i < 3
%                 title('Knee','Fontsize',label_fontsize)
%             end
%             if (i > 2 && i < 5)
%                 title('Hip flexion','Fontsize',label_fontsize)
%             end
%             if (i > 4 && i < 7)
%                 title('Hip adduction','Fontsize',label_fontsize)
%             end
%             if i == 1
%                 ylabel('Max R angle','Fontsize',label_fontsize);
%             end
%             if i == 13
%                 ylabel('Max R torque','Fontsize',label_fontsize);
%             end
%             if i == 7 || i == 19
%                 ylabel('R shift (%)','Fontsize',label_fontsize);
%             end
%             if i == 19 || i == 21 || i == 23
%                 xlabel('k tmt (Nm/rad)','Fontsize',label_fontsize);
%             end
%             if i == 20 || i == 22 || i == 24
%                 xlabel('cWL (/°)','Fontsize',label_fontsize);
%             end
%             
%         end
% 
%     end
%     
%     if boolFirst
%         lh=legend('-DynamicLegend','location','east');
%         lh.Interpreter = 'none';
%         lhPos = lh.Position;
%         lhPos(1) = lhPos(1)+0.2;
%         set(lh,'position',lhPos);
%     end
%     
else
    
end

end