function [] = PlotResults_3DSim_tmt(ResultsFile,Cs,LegName,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


% Notes: We have to update the computations of the biological joint moments
% because this is currently based on subtracting the exo torque from the
% total ankle moment. Since we actuate more than 1 joint, we will have to
% do this computation based on R.Exodiff_id

% varargin:
%   (1) handle figure
%   (2) value x-axis (COT graph)
%   (3) label x-axis (COT graph)



set(0,'defaultTextInterpreter','none');

if ~exist('LegName','var')
    [~,filename,~]= fileparts(ResultsFile);
    LegName = filename;
end



if exist(ResultsFile,'file')
    load(ResultsFile);
    
    % option to show only simulation results, no measurement data
    if strcmp(varargin{length(varargin)},'no_meas_data') || strcmp(varargin{length(varargin)},'none')
        md = 0;
        vv = 0;
    else
        md = 1;
        if strcmp(varargin{length(varargin)},'act')
            type = 'Active';
            vv = 0;
        elseif strcmp(varargin{length(varargin)},'pas')
            type = 'Passive';
            vv = 0;
        elseif strcmp(varargin{length(varargin)},'norm')
            type = 'Normal';
            vv = 0;
        else
            type = 'Normal';
        end
    end
    
    if length(varargin)<2 || (~vv && length(varargin)<3)
        xParam = R.Sopt.ExoScale;
        xParamLab = 'Level assistance';
    else
        xParam = varargin{2};
        xParamLab = varargin{3};
    end
    
    
    
    % determine if the model has tmt joint
    has_no_tmt = ~isfield(R.S,'tmt') || isempty(R.S.tmt) || R.S.tmt == 0;
    
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
            tab6 = h.Parent.Children(1).Children(1).Children(6);
            tab7 = h.Parent.Children(1).Children(1).Children(7);
            tab8 = h.Parent.Children(1).Children(1).Children(8);
            tab9 = h.Parent.Children(1).Children(1).Children(9);
            boolFirst = 0;
        else
            h = varargin{1};
            hTabGroup = uitabgroup;
            tab1 = uitab(hTabGroup, 'Title', 'Kinematics');
            tab2 = uitab(hTabGroup, 'Title', 'Kinetics');
            tab3 = uitab(hTabGroup, 'Title', 'COT');
            tab4 = uitab(hTabGroup, 'Title', 'ExoInfo');
            tab5 = uitab(hTabGroup, 'Title', 'CalfM');
            tab6 = uitab(hTabGroup, 'Title', 'Ground reaction force');
            tab7 = uitab(hTabGroup, 'Title', 'Objective Function');
            tab8 = uitab(hTabGroup, 'Title', 'Ankle detailed');
            tab9 = uitab(hTabGroup, 'Title', 'SpatioTemporal');
            h.Name = 'Sim3D_Results';
            set(h,'Color','w');
        end
    else
        h = figure();
        h.Name = 'Sim3D_Results';
        hTabGroup = uitabgroup;
        tab1 = uitab(hTabGroup, 'Title', 'Kinematics');
        tab2 = uitab(hTabGroup, 'Title', 'Kinetics');
        tab3 = uitab(hTabGroup, 'Title', 'COT');
        tab4 = uitab(hTabGroup, 'Title', 'ExoInfo');
        tab5 = uitab(hTabGroup, 'Title', 'CalfM');
        tab6 = uitab(hTabGroup, 'Title', 'Ground reaction force');
        tab7 = uitab(hTabGroup, 'Title', 'Objective Function');
        tab8 = uitab(hTabGroup, 'Title', 'Ankle detailed');
        tab9 = uitab(hTabGroup, 'Title', 'SpatioTemporal');
        set(h,'Color','w');
    end
    
    %% Helper names
%     joints_ref = {'pelvis_tilt','pelvis_list','pelvis_rotation',...
%         'hip_flexion','hip_adduction','hip_rotation',...
%         'knee_angle','ankle_angle','subtalar_angle','tmt_angle','mtp_angle',...
%         'lumbar_extension','lumbar_bending','lumbar_rotation',...
%         'arm_flex','arm_add','arm_rot','elbow_flex'};
    joints_ref = {'pelvis_tilt','pelvis_list','pelvis_rotation',...
        'hip_flexion_r','hip_adduction_r','hip_rotation_r',...
        'knee_angle_r','ankle_angle_r','subtalar_angle_r','tmt_angle_r','mtp_angle_r',...
        'lumbar_extension','lumbar_bending','lumbar_rotation',...
        'arm_flex_r','arm_add_r','arm_rot_r','elbow_flex_r'};
    joints_tit = {'Pelvis tilt','Pelvis list','Pelvis rotation','Pelvis tx',...
        'Pelvis ty','Pelvis tz','Hip flexion L','Hip adduction L',...
        'Hip rotation L','Hip flexion R','Hip adduction R','Hip rotation R',...
        'Knee L','Knee R','Ankle L','Ankle R',...
        'Subtalar L','Subtalar R','TMT L','TMT R','MTP L','MTP R',...
        'Lumbar extension','Lumbar bending','Lumbar rotation',...
        'Arm flexion L','Arm adduction L','Arm rotation L',...
        'Arm flexion R','Arm adduction R','Arm rotation R',...
        'Elbow flexion L','Elbow flexion R'};
    
    %% Plot default figure kinematics (Antoine his code)
    axes('parent', tab1);
    
    if boolFirst && md
        subject = R.S.subject;
%         [pathHere,~,~] = fileparts(mfilename('fullpath'));
%         [pathrepo,~,~] = fileparts(pathHere);
% 
%         pathReferenceData = [pathrepo,'/ExperimentalData'];
%         load([pathReferenceData,'/ExperimentalData.mat'],'ExperimentalData');
%         Qref = ExperimentalData.Q.(subject).Qs;

        
        % load data Pog_s1 from struct saved during ...\Analyze_ExoData\Batch\BatchScript_LatexReport.m
        load('D:\school\WTK\thesis\model\3dpredictsim\Data\Pog_s1.mat','Dat');
        Qref = Dat.(type).gc;

        
    end
    
    if has_no_tmt
        idx_Qs = [1,2,3,10,11,12,14,16,18,20,21,22,23,27,28,29,31];
    else
        idx_Qs = [1,2,3,10,11,12,14,16,18,20,22,23,24,25,29,30,31,33];
    end
    idx_title = [1,2,3,10,11,12,14,16,18,20,22,23,24,25,29,30,31,33];
    
    j = 0;
    label_fontsize  = 12;
    line_linewidth  = 2;
    for i = 1:length(idx_title)
        subplot(3,6,i)
        x = 1:(100-1)/(size(R.Qs,1)-1):100;
        % Experimental data
        if  boolFirst == 1 && md
            idx_jref = strcmp(Qref.colheaders,joints_ref{i});
            if sum(idx_jref) == 1
%                 meanPlusSTD = Qref.(subject).Qs.mean(:,idx_jref) + 2*Qref.(subject).Qs.std(:,idx_jref);
%                 meanMinusSTD = Qref.(subject).Qs.mean(:,idx_jref) - 2*Qref.(subject).Qs.std(:,idx_jref);
                meanPlusSTD = (Qref.Qall_mean(:,idx_jref) + 2*Qref.Qall_std(:,idx_jref)).*180/pi;
                meanMinusSTD = (Qref.Qall_mean(:,idx_jref) - 2*Qref.Qall_std(:,idx_jref)).*180/pi;
                
                stepQ = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
                intervalQ = 1:stepQ:size(R.Qs,1);
                sampleQ = 1:size(R.Qs,1);
                meanPlusSTD = interp1(intervalQ,meanPlusSTD,sampleQ);
                meanMinusSTD = interp1(intervalQ,meanMinusSTD,sampleQ);

                hold on
                fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k','DisplayName',['MoCap ' subject '(' type ')']);
                alpha(.25);
            end
        end
        
        % Simulation results
        x = 1:(100-1)/(size(R.Qs,1)-1):100;
        hold on;
        if has_no_tmt && (strcmp(joints_ref{i},'tmt_angle') || strcmp(joints_ref{i},'tmt_angle_r'))
            % skip this plot
        else
            j=j+1;
        
            if i == length(idx_title)
                 plot(x,R.Qs(:,idx_Qs(j)),'color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
            else
                plot(x,R.Qs(:,idx_Qs(j)),'color',Cs,'linewidth',line_linewidth);
            end
        end
        
        % Plot settings
        if boolFirst
            set(gca,'Fontsize',label_fontsize);
            title(joints_tit{idx_title(i)},'Fontsize',label_fontsize);
            % Y-axis
            if i == 1 || i == 7 || i == 13 
                ylabel('Angle (°)','Fontsize',label_fontsize);
            end
            % X-axis
            L = get(gca,'XLim');
            NumTicks = 2;
            if i > 12
                set(gca,'XTick',linspace(L(1),L(2),NumTicks))
                xlabel('Gait cycle (%)','Fontsize',label_fontsize);
            else
                set(gca,'XTick',[]);
            end
        end
    end
    
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
        lhPos = lh.Position;
%         lhPos(1) = lhPos(1)+0.2;
        set(lh,'position',lhPos);
    end
    
    %% Plot joint kinetics
    axes('parent', tab2);
    j = 0;
    label_fontsize  = 12;
    line_linewidth  = 2;
    for i = 1:length(idx_title)
        subplot(3,6,i)
        x = 1:(100-1)/(size(R.Qs,1)-1):100;
        % Experimental data
        if  boolFirst == 1 && md
%             IDref = ExperimentalData.Torques;
%             idx_jref = strcmp(IDref.(subject).colheaders,joints_ref{i});
            idx_jref = strcmp(Qref.colheaders,joints_ref{i});
            if sum(idx_jref) == 1
%                 meanPlusSTD = IDref.(subject).mean(:,idx_jref) + 2*IDref.(subject).std(:,idx_jref);
%                 meanMinusSTD = IDref.(subject).mean(:,idx_jref) - 2*IDref.(subject).std(:,idx_jref);
                meanPlusSTD = Qref.Tall_bio_mean(:,idx_jref) + 2*Qref.Tall_bio_std(:,idx_jref);
                meanMinusSTD = Qref.Tall_bio_mean(:,idx_jref) - 2*Qref.Tall_bio_std(:,idx_jref);
                
                stepID = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
                intervalID = 1:stepID:size(R.Qs,1);
                sampleID = 1:size(R.Qs,1);
                meanPlusSTD = interp1(intervalID,meanPlusSTD,sampleID);
                meanMinusSTD = interp1(intervalID,meanMinusSTD,sampleID);
                hold on
                fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k','DisplayName',['MoCap ' subject]);
                alpha(.25);
            end
        end
        
        % Simulation results
        x = 1:(100-1)/(size(R.Qs,1)-1):100;
        hold on;
        if has_no_tmt && (strcmp(joints_ref{i},'tmt_angle') || strcmp(joints_ref{i},'tmt_angle_r'))
            % skip this plot
        else
            j=j+1;
        
            if i == length(idx_title)
                plot(x,R.Tid(:,idx_Qs(j)),'color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
            else
                plot(x,R.Tid(:,idx_Qs(j)),'color',Cs,'linewidth',line_linewidth);
            end
        end
        % Plot settings
        if boolFirst
            % Plot settings
            set(gca,'Fontsize',label_fontsize);
            title(joints_tit{idx_title(i)},'Fontsize',label_fontsize);
            % Y-axis
            if i == 1 || i == 7 ||i == 13
                ylabel('Torque (Nm)','Fontsize',label_fontsize);
            end
            % X-axis
            L = get(gca,'XLim');
            NumTicks = 2;
            if i > 12
                set(gca,'XTick',linspace(L(1),L(2),NumTicks))
                xlabel('Gait cycle (%)','Fontsize',label_fontsize);
            else
                set(gca,'XTick',[]);
            end
        end
    end
    
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
        lhPos = lh.Position;
%         lhPos(1) = lhPos(1)+0.2;
        set(lh,'position',lhPos);
    end
    
    
    %% Plot general information
    axes('parent', tab3);
    
    % Plot COT as a function of exo assitance
    subplot(2,2,1); hold on;
    % experimental data
    darkgray = [0.5 0.5 0.5];
    lightgray = [0.1 0.1 0.1];
    if boolFirst && md
        plot(0,Dat.Passive.COT,'o','Color',darkgray,'MarkerFaceColor','k','DisplayName','Data exo passive');
        plot(0,Dat.Normal.COT,'o','Color',darkgray,'MarkerFaceColor',darkgray,'DisplayName','Data normal shoes');
        plot(1,Dat.Active.COT,'o','Color',lightgray,'MarkerFaceColor',lightgray,'DisplayName','Data exo active');
    end
    % simulation results
    plot(xParam,R.COT,'o','Color',Cs,'MarkerFaceColor',Cs,'DisplayName',LegName);
    ylabel('COT');  xlabel(xParamLab);
    title('Cost of Transport');
    xlim([-0.5,1.5])
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    
    % Plot stride frequency
    subplot(2,2,2); hold on;
    if boolFirst && md
        plot(0,Dat.Normal.Step.StrideFreq_mean,'o','Color',darkgray,'MarkerFaceColor',darkgray);
        errorbar(0,Dat.Normal.Step.StrideFreq_mean,Dat.Normal.Step.StrideFreq_std,'Color',darkgray,'LineWidth',1.5);
        plot(1,Dat.Active.Step.StrideFreq_mean,'o','Color',lightgray,'MarkerFaceColor',lightgray);
        errorbar(1,Dat.Active.Step.StrideFreq_mean,Dat.Normal.Step.StrideFreq_std,'Color',lightgray,'LineWidth',1.5);
    end
    dt = R.t(end);
    plot(xParam,1./dt,'o','Color',Cs,'MarkerFaceColor',Cs);
    ylabel('Stride Frequency'); xlabel(xParamLab);
    title('Stride Frequency');
    xlim([-0.5,1.5])
    
    subplot(2,2,3);  hold on;
    plot(xParam,R.dt_exoShift,'o','Color',Cs,'MarkerFaceColor',Cs);
    ylabel('Time shift exo torque[s]');  xlabel(xParamLab);
    
    
    subplot(2,2,4); hold on;
    plot(R.T_exo(:,2),'-','Color',Cs);
    ylabel('Exo Moment [Nm]');  xlabel('% stride');
    title('Right');
    
    %% Plot Torque information
    % update this here
    boolActuation = 0;
    if isfield(R,'Exodiff_id')
        boolActuation = 1;
    end
    axes('parent', tab4);
    iAnkle = strcmp(R.colheaders.joints,'ankle_angle_r');
    iSubtalar = strcmp(R.colheaders.joints,'subtalar_angle_r');
    subplot(3,2,1)
    if boolActuation
        plot(R.Exodiff_id(:,iAnkle),'-','Color',Cs); hold on;
    else
        plot(R.T_exo(:,2),'-','Color',Cs); hold on;
    end
    ylabel('Exo Moment - Ankle [Nm]');  xlabel('% stride');
    
    subplot(3,2,2);
    if boolActuation
        plot(R.Exodiff_id(:,iSubtalar),'-','Color',Cs); hold on;
    else
    end
    ylabel('Exo Moment- Subtalar [Nm]'); xlabel('% stride');
    
    subplot(3,2,3)
    plot(R.Tid(:,iAnkle),'-','Color',Cs); hold on;
    ylabel('Ankle moment [Nm]'); xlabel('% stride');
    
    subplot(3,2,4)
    plot(R.Tid(:,iSubtalar),'-','Color',Cs); hold on;
    ylabel('Subtalar moment [Nm]'); xlabel('% stride');
    
    subplot(3,2,5)
    if boolActuation
        plot(R.Tid(:,iAnkle)-R.Exodiff_id(:,iAnkle),'-','Color',Cs); hold on;
    else
        plot(R.Tid(:,iAnkle)-R.T_exo(:,2),'-','Color',Cs); hold on;
    end
    ylabel('Biological ankle moment [Nm]'); xlabel('% stride');
    title('Left');
    
    subplot(3,2,6)
    if boolActuation
        plot(R.Tid(:,iSubtalar)-R.Exodiff_id(:,iSubtalar),'-','Color',Cs,'DisplayName',LegName); hold on;
    else
        plot(R.Tid(:,iSubtalar),'-','Color',Cs,'DisplayName',LegName); hold on;
    end
    ylabel('Biological subtalar moment [Nm]'); xlabel('% stride');
    title('Left');
    
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    
    %% Plot ankle muscle energetics
    axes('parent', tab5);
    iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
    iGas = find(strcmp(R.colheaders.muscles,'lat_gas_r'));
 
    subplot(5,2,1); hold on;
    plot(R.a(:,iSol),'-','Color',Cs,'DisplayName',LegName);  title('Soleus');
    xlabel('% stride'); ylabel('activity');
    
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    
    subplot(5,2,2); hold on;
%     if boolFirst && md
%         plot(gas_act_normal,'-','Color',darkgray);
%         plot(gas_act_active,'-','Color',lightgray);
%     end
    plot(R.a(:,iGas),'-','Color',Cs); title('Gastrocnemius');
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
    plot(R.FT(:,iGas),'-','Color',Cs,'DisplayName',LegName); hold on;
    xlabel('% stride'); ylabel('Norm muscle force');
    
    
    
    %% Ground reaction force
    axes('parent', tab6);
    for i=1:6
        subplot(2,3,i);
        l = plot(R.GRFs(:,i),'-','Color',Cs); hold on;
        title(R.colheaders.GRF{i});
        xlabel('% stride');
    end
    l.DisplayName = LegName;
     if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    
    %% Objective function
    axes('parent', tab7);
    if isfield(R,'Obj')
        Fields = fieldnames(R.Obj);
        nf = length(Fields);
        for i= 1:nf
            subplot(3,4,i)
            l = plot(xParam,R.Obj.(Fields{i}),'o','Color',Cs,'MarkerFaceColor',Cs); hold on;
            xlabel(xParamLab);
            title(Fields{i});
        end
        l.DisplayName = LegName;
        if boolFirst
            lh=legend('-DynamicLegend','location','east');
            lh.Interpreter = 'none';
        end
    end
    
    
    %% Plot ankle muscle energetics
    axes('parent', tab8);
    iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
    iGas = find(strcmp(R.colheaders.muscles,'lat_gas_r'));
    iGas2 = find(strcmp(R.colheaders.muscles,'med_gas_r'));
    iTib = find(strcmp(R.colheaders.muscles,'tib_ant_r'));
    if isempty(iGas)
        iGas = find(strcmp(R.colheaders.muscles,'gaslat_r'));
    end
    if isempty(iGas2)
        iGas2 = find(strcmp(R.colheaders.muscles,'gasmed_r'));
    end
    if isempty(iTib)
        iTib = find(strcmp(R.colheaders.muscles,'tibant_r'));
    end
    mVect = {'Soleus','Gas-lat','Gas-med','Tib-ant'};
    
    iM = [iSol iGas iGas2 iTib];
    
    for i=1:4
        subplot(5,4,i)
        plot(R.a(:,iM(i)),'-','Color',Cs); hold on; title(mVect{i});
        xlabel('% stride'); ylabel('activity');
        
        subplot(5,4,i+4)
        plot(R.MetabB.Etot(:,iM(i)),'-','Color',Cs); hold on;
        xlabel('% stride'); ylabel('Muscle metab power');
        
        subplot(5,4,i+8)
        plot(R.lMtilde(:,iM(i)),'-','Color',Cs); hold on;
        xlabel('% stride'); ylabel('Norm fiber length');
        
        subplot(5,4,i+12)
        plot(R.FT(:,iM(i)),'-','Color',Cs); hold on;
        xlabel('% stride'); ylabel('Norm muscle force');
        
    end
    % plot (biological) joint moments
    subplot(5,4,17)
    plot(R.Tid(:,strcmp(R.colheaders.joints,'ankle_angle_r')),'-','Color',Cs); hold on;
    ylabel('Ankle moment [Nm]'); xlabel('% stride');
    
    subplot(5,4,18)
    plot(R.Tid(:,strcmp(R.colheaders.joints,'ankle_angle_r'))-R.T_exo(:,2),'-','Color',Cs); hold on;
    ylabel('Muscle ankle [Nm]'); xlabel('% stride');
    
    subplot(5,4,19)
    plot(R.Tid(:,strcmp(R.colheaders.joints,'subtalar_angle_r')),'-','Color',Cs); hold on;
    ylabel('subtalar moment [Nm]'); xlabel('% stride');
    
    subplot(5,4,20)
    l = plot(R.Tid(:,strcmp(R.colheaders.joints,'mtp_angle_r')),'-','Color',Cs); hold on;
    ylabel('mtp moment [Nm]'); xlabel('% stride');
    l.DisplayName = LegName;
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    
    %% Spatiotemporal results
    
    axes('parent', tab9);
    
    subplot(2,3,1); hold on;
    if boolFirst && md
        plot(0,Dat.Passive.speed./Dat.Normal.Step.StrideFreq_mean,'o','Color','k','MarkerFaceColor','k','DisplayName','Data normal shoes');
        plot(0,Dat.Normal.speed./Dat.Normal.Step.StrideFreq_mean,'o','Color',darkgray,'MarkerFaceColor',darkgray,'DisplayName','Data normal shoes');
        plot(1,Dat.Active.speed./Dat.Active.Step.StrideFreq_mean,'o','Color',lightgray,'MarkerFaceColor',lightgray,'DisplayName','Data exo active');
    end
    plot(xParam,R.StrideLength,'o','Color',Cs,'MarkerFaceColor',Cs)
    ylabel('Stride length [m]');  xlabel(xParamLab);
    xlim([-0.5,1.5])
    
    subplot(2,3,2); hold on;
    if boolFirst && md
        mean_norm = nanmean(Dat.Normal.StepWidth);
        mean_act = nanmean(Dat.Active.StepWidth);
        mean_pas = nanmean(Dat.Passive.StepWidth);
%         plot(0,mean_norm,'o','Color',darkgray,'MarkerFaceColor',darkgray);
%         plot(1,mean_act,'o','Color',lightgray,'MarkerFaceColor',lightgray);
        errorbar(0,mean_pas,nanstd(Dat.Passive.StepWidth),'Color','k','LineWidth',1.5);
        errorbar(0,mean_norm,nanstd(Dat.Normal.StepWidth),'Color',darkgray,'LineWidth',1.5);
        errorbar(1,mean_act,nanstd(Dat.Active.StepWidth),'Color',lightgray,'LineWidth',1.5);
    end
    if isfield(R,'StepWidth_COP')
        plot(xParam,R.StepWidth_COP,'o','Color',Cs,'MarkerFaceColor',Cs); 
    end
    ylabel('Stride width [m]');  xlabel(xParamLab);
    xlim([-0.5,1.5])
    
    subplot(2,3,3); hold on;
    if boolFirst && md
%         plot(0,Dat.Normal.Step.StrideFreq_mean,'o','Color',darkgray,'MarkerFaceColor',darkgray);
%         plot(1,Dat.Active.Step.StrideFreq_mean,'o','Color',lightgray,'MarkerFaceColor',lightgray);
        errorbar(0,Dat.Passive.Step.StrideFreq_mean,Dat.Passive.Step.StrideFreq_std,'Color','k','LineWidth',1.5);
        errorbar(0,Dat.Normal.Step.StrideFreq_mean,Dat.Normal.Step.StrideFreq_std,'Color',darkgray,'LineWidth',1.5);
        errorbar(1,Dat.Active.Step.StrideFreq_mean,Dat.Normal.Step.StrideFreq_std,'Color',lightgray,'LineWidth',1.5);
    end    
    plot(xParam,1./(R.tf_step*2),'o','Color',Cs,'MarkerFaceColor',Cs);
    ylabel('stride frequency [Hz]');  xlabel(xParamLab);
    xlim([-0.5,1.5])
    
    % probably mistake in std calculation data
    subplot(2,3,4); hold on;
    if boolFirst && md
%         plot(0,Dat.Normal.Step.PercStance,'o','Color',darkgray,'MarkerFaceColor',darkgray,'DisplayName','Data normal shoes');
%         plot(1,Dat.Active.Step.PercStance,'o','Color',lightgray,'MarkerFaceColor',lightgray,'DisplayName','Data exo active');
        errorbar(0,Dat.Passive.Step.PercStance,mean(nanstd(Dat.Passive.Step.PercStance_n)),'Color','k','LineWidth',1.5);
        errorbar(0,Dat.Normal.Step.PercStance,mean(nanstd(Dat.Normal.Step.PercStance_n)),'Color',darkgray,'LineWidth',1.5);
        errorbar(1,Dat.Active.Step.PercStance,mean(nanstd(Dat.Active.Step.PercStance_n)),'Color',lightgray,'LineWidth',1.5);
    end
    plot(xParam,R.Event.Stance,'o','Color',Cs,'MarkerFaceColor',Cs);
    ylabel('% stance');  xlabel(xParamLab);
    xlim([-0.5,1.5])
    
    subplot(2,3,5); hold on;
    if boolFirst && md
%         plot(0,Dat.Normal.Step.PercSwing,'o','Color',darkgray,'MarkerFaceColor',darkgray,'DisplayName','Data normal shoes');
%         plot(1,Dat.Active.Step.PercSwing,'o','Color',lightgray,'MarkerFaceColor',lightgray,'DisplayName','Data exo active');
        errorbar(0,Dat.Passive.Step.PercSwing,mean(nanstd(Dat.Passive.Step.PercSwing_n)),'Color','k','LineWidth',1.5,'DisplayName','Data exo passive');
        errorbar(0,Dat.Normal.Step.PercSwing,mean(nanstd(Dat.Normal.Step.PercSwing_n)),'Color',darkgray,'LineWidth',1.5,'DisplayName','Data normal shoes');
        errorbar(1,Dat.Active.Step.PercSwing,mean(nanstd(Dat.Active.Step.PercSwing_n)),'Color',lightgray,'LineWidth',1.5,'DisplayName','Data exo active');
    end
    plot(xParam,R.Event.Swing,'o','Color',Cs,'MarkerFaceColor',Cs);
    ylabel('% swing');  xlabel(xParamLab);
    xlim([-0.5,1.5])
    
    subplot(2,3,6); hold on;
    if boolFirst && md
%         plot(0,Dat.Normal.Step.PercDS,'o','Color',darkgray,'MarkerFaceColor',darkgray,'DisplayName','Data normal shoes');
%         plot(1,Dat.Active.Step.PercDS,'o','Color',lightgray,'MarkerFaceColor',lightgray,'DisplayName','Data exo active');
        errorbar(0,Dat.Passive.Step.PercDS,mean(nanstd(Dat.Passive.Step.PercDS_n)),'Color','k','LineWidth',1.5,'DisplayName','Data exo passive');
        errorbar(0,Dat.Normal.Step.PercDS,mean(nanstd(Dat.Normal.Step.PercDS_n)),'Color',darkgray,'LineWidth',1.5,'DisplayName','Data normal shoes');
        errorbar(1,Dat.Active.Step.PercDS,mean(nanstd(Dat.Active.Step.PercDS_n)),'Color',lightgray,'LineWidth',1.5,'DisplayName','Data exo active');
    end
    l = plot(xParam,R.Event.DS,'o','Color',Cs,'MarkerFaceColor',Cs);
    ylabel('% double support');  xlabel(xParamLab);
    xlim([-0.5,1.5])
    
    l.DisplayName = LegName;
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    
else
    warning(['File not found: ' ResultsFile]);
end

end

