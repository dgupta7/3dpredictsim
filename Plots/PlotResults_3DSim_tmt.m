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
    subject = R.S.subject;
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
        elseif strcmp(varargin{length(varargin)},'Fal_s1')
            type = 'Normal';
            vv = 0;
            subject = 'Fal_s1';
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
    
    if strcmp(R.S.subject,'subject1')
        subject = 'Fal_s1';
    end
    
    % determine if the model has tmt joint
    has_no_tmt = ~isfield(R.S,'tmt') || isempty(R.S.tmt) || R.S.tmt == 0;
    has_no_mtj = ~isfield(R.S,'mtj') || isempty(R.S.mtj) || R.S.mtj == 0;
    
    
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
            tab10 = h.Parent.Children(1).Children(1).Children(10);
            tab11 = h.Parent.Children(1).Children(1).Children(11);
            tab12 = h.Parent.Children(1).Children(1).Children(12);
            tab13 = h.Parent.Children(1).Children(1).Children(13);
            tab14 = h.Parent.Children(1).Children(1).Children(14);
            tab15 = h.Parent.Children(1).Children(1).Children(15);
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
            tab7 = uitab(hTabGroup, 'Title', 'GRF detailed');
            tab8 = uitab(hTabGroup, 'Title', 'Objective Function');
            tab9 = uitab(hTabGroup, 'Title', 'Ankle detailed');
            tab10 = uitab(hTabGroup, 'Title', 'Mtp detailed');
            tab11 = uitab(hTabGroup, 'Title', 'SpatioTemporal');
            tab12 = uitab(hTabGroup, 'Title', 'Windlass');
            tab13 = uitab(hTabGroup, 'Title', 'Foot power');
            tab14 = uitab(hTabGroup, 'Title', 'Foot work');
            tab15 = uitab(hTabGroup, 'Title', 'Exo assistance');
            
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
        tab7 = uitab(hTabGroup, 'Title', 'GRF detailed');
        tab8 = uitab(hTabGroup, 'Title', 'Objective Function');
        tab9 = uitab(hTabGroup, 'Title', 'Ankle detailed');
        tab10 = uitab(hTabGroup, 'Title', 'Mtp detailed');
        tab11 = uitab(hTabGroup, 'Title', 'SpatioTemporal');
        tab12 = uitab(hTabGroup, 'Title', 'Windlass');
        tab13 = uitab(hTabGroup, 'Title', 'Foot power');
        tab14 = uitab(hTabGroup, 'Title', 'Foot work');
        tab15 = uitab(hTabGroup, 'Title', 'Exo assistance');
        
        set(h,'Color','w');
    end
    
    %% Helper names
    if strcmp(subject,'Fal_s1')
        joints_ref = {'pelvis_tilt','pelvis_list','pelvis_rotation',...
            'hip_flexion','hip_adduction','hip_rotation',...
            'knee_angle','ankle_angle','subtalar_angle','mtj_angle','tmt_angle','mtp_angle',...
            'lumbar_extension','lumbar_bending','lumbar_rotation',...
            'arm_flex','arm_add','arm_rot','elbow_flex'};
    else
        joints_ref = {'pelvis_tilt','pelvis_list','pelvis_rotation',...
            'hip_flexion_r','hip_adduction_r','hip_rotation_r',...
            'knee_angle_r','ankle_angle_r','subtalar_angle_r','mtj_angle_r','tmt_angle_r','mtp_angle_r',...
            'lumbar_extension','lumbar_bending','lumbar_rotation',...
            'arm_flex_r','arm_add_r','arm_rot_r','elbow_flex_r'};
    end
    joints_tit = {'Pelvis tilt','Pelvis list','Pelvis rotation','Pelvis tx',...
        'Pelvis ty','Pelvis tz','Hip flexion L','Hip adduction L',...
        'Hip rotation L','Hip flexion R','Hip adduction R','Hip rotation R',...
        'Knee L','Knee R','Ankle L','Ankle R','Subtalar L','Subtalar R',...
        'Midtarsal L','Midtarsal R','Tarsometatarsal L','Tarsometatarsal R','MTP L','MTP R',...
        'Lumbar extension','Lumbar bending','Lumbar rotation',...
        'Arm flexion L','Arm adduction L','Arm rotation L',...
        'Arm flexion R','Arm adduction R','Arm rotation R',...
        'Elbow flexion L','Elbow flexion R'};
    
    %% Plot default figure kinematics (Antoine his code)
    axes('parent', tab1);
    
    if boolFirst && md
        [pathHere,~,~] = fileparts(mfilename('fullpath'));
        [pathRepo,~,~] = fileparts(pathHere);
        try
            if strcmp(subject,'Fal_s1')
                load([pathRepo '\Data\Fal_s1.mat'],'Dat');
            else
                load([pathRepo '\Data\Pog_s1.mat'],'Dat');
            end
            Qref = Dat.(type).gc;
        catch
            md = 0;
            disp('Reference data not found')

        end
        
    end
    
    if has_no_tmt && has_no_mtj
        idx_Qs = [1,2,3,10,11,12,14,16,18,20,21,22,23,27,28,29,31];
    else
        idx_Qs = [1,2,3,10,11,12,14,16,18,20,22,23,24,25,29,30,31,33];
    end
    idx_title = [1,2,3,10,11,12,14,16,18,20,22,24,25,26,27,31,32,33,35];
    
    j = 0;
    label_fontsize  = 12;
    line_linewidth  = 0.5;
    for i = 1:length(idx_title)
        subplot(3,7,i)
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
        if (has_no_tmt && strcmp(joints_tit{idx_title(i)},'Tarsometatarsal R')) || ...
                (has_no_mtj && strcmp(joints_tit{idx_title(i)},'Midtarsal R'))
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
            if i == 1 || i == 8 || i == 15 
                ylabel('Angle (�)','Fontsize',label_fontsize);
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
        lh=legend('-DynamicLegend','location','west');
        lh.Interpreter = 'none';
        lhPos = lh.Position;
        lhPos(1) = lhPos(1)+0.1;
        set(lh,'position',lhPos);
    end
    
    %% Plot joint kinetics
    axes('parent', tab2);
    j = 0;
    label_fontsize  = 12;
    line_linewidth  = 0.5;
    for i = 1:length(idx_title)
        subplot(3,7,i)
        x = 1:(100-1)/(size(R.Qs,1)-1):100;
        % Experimental data
        if  boolFirst == 1 && md
%             IDref = ExperimentalData.Torques;
%             idx_jref = strcmp(IDref.(subject).colheaders,joints_ref{i});
            idx_jref = strcmp(Qref.colheaders,joints_ref{i});
            if sum(idx_jref) == 1
%                 meanPlusSTD = IDref.(subject).mean(:,idx_jref) + 2*IDref.(subject).std(:,idx_jref);
%                 meanMinusSTD = IDref.(subject).mean(:,idx_jref) - 2*IDref.(subject).std(:,idx_jref);
                meanPlusSTD = Qref.Tall_mean(:,idx_jref) + 2*Qref.Tall_std(:,idx_jref);
                meanMinusSTD = Qref.Tall_mean(:,idx_jref) - 2*Qref.Tall_std(:,idx_jref);
                
                stepID = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
                intervalID = 1:stepID:size(R.Qs,1);
                sampleID = 1:size(R.Qs,1);
                meanPlusSTD = interp1(intervalID,meanPlusSTD,sampleID);
                meanMinusSTD = interp1(intervalID,meanMinusSTD,sampleID);
                hold on
                fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k','DisplayName',['MoCap ' subject '(' type ')']);
                alpha(.25);
            end
        end
        
        % Simulation results
        x = 1:(100-1)/(size(R.Qs,1)-1):100;
        hold on;
        if (has_no_tmt && strcmp(joints_tit{idx_title(i)},'Tarsometatarsal R')) || ...
                (has_no_mtj && strcmp(joints_tit{idx_title(i)},'Midtarsal R'))
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
            if i == 1 || i == 8 ||i == 15
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
        lh=legend('-DynamicLegend','location','west');
        lh.Interpreter = 'none';
        lhPos = lh.Position;
        lhPos(1) = lhPos(1)+0.1;
        set(lh,'position',lhPos);
    end
    
    
    %% Plot general information
    axes('parent', tab3);
    
    % Plot COT as a function of exo assitance
    subplot(2,2,1); hold on;
    % experimental data
    darkgray = [0.5 0.5 0.5];
    lightgray = [0.1 0.1 0.1];
    if boolFirst && md && ~strcmp(subject,'Fal_s1') 
        plot(0,Dat.Normal.COT,'o','Color',darkgray,'MarkerFaceColor',darkgray,'DisplayName','Data normal shoes');
        plot(0,Dat.Passive.COT,'o','Color',darkgray,'MarkerFaceColor','k','DisplayName','Data exo passive');
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
    if boolFirst && md && ~strcmp(subject,'Fal_s1')
        plot(0,Dat.Normal.Step.StrideFreq_mean,'o','Color',darkgray,'MarkerFaceColor',darkgray);
        errorbar(0,Dat.Normal.Step.StrideFreq_mean,Dat.Normal.Step.StrideFreq_std,'Color',darkgray,'LineWidth',1.5);
        plot(1,Dat.Active.Step.StrideFreq_mean,'o','Color',lightgray,'MarkerFaceColor',lightgray);
        errorbar(1,Dat.Active.Step.StrideFreq_mean,Dat.Normal.Step.StrideFreq_std,'Color',lightgray,'LineWidth',1.5);
    end
%     dt = R.t(end);
%     plot(xParam,1./dt,'o','Color',Cs,'MarkerFaceColor',Cs);
    plot(xParam,1./(R.tf_step*2),'o','Color',Cs,'MarkerFaceColor',Cs);
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
    iTmt = strcmp(R.colheaders.joints,'tmt_angle_r');
    
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
    
    subplot(5,2,2); hold on;
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
    
    if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    
    %% Ground reaction force
    axes('parent', tab6);
    
    if boolFirst && md 
        for i=1:3
            subplot(2,3,i)
            hold on
            if strcmp(subject,'Fal_s1')
                plot(Dat.(type).gc.GRF.Fmean(:,i),'-k');
            else
                plot(Dat.(type).gc.GRF.Fmean(:,i)/(R.body_mass*9.81)*100,'-k');
            end
        end
    end
    for i=1:6
        subplot(2,3,i)
        hold on
        l = plot(R.GRFs(:,i),'-','Color',Cs);
        title(R.colheaders.GRF{i});
        xlabel('% stride');
        ylabel('% body weight')
    end
    l.DisplayName = LegName;
     if boolFirst
        lh=legend('-DynamicLegend','location','east');
        lh.Interpreter = 'none';
    end
    
%% GRF detailed
    axes('parent', tab7);
    
    if boolFirst && md
        for i=1:3
            subplot(3,4,i)
            hold on
            if strcmp(subject,'Fal_s1')
                plot(Dat.(type).gc.GRF.Fmean(:,i),'-k');
            else
                plot(Dat.(type).gc.GRF.Fmean(:,i)/(R.body_mass*9.81)*100,'-k');
            end
        end
    end
    for i=1:3
        subplot(3,4,i)
        hold on
        l = plot(R.GRFs(:,i),'-','Color',Cs);
        title(R.colheaders.GRF{i});
        xlabel('% stride');
        ylabel('% body weight')
    end
    l.DisplayName = LegName;
    

    if isfield(R,'GRFs_separate') && ~isempty(R.GRFs_separate)
        for i=1:3
            subplot(3,4,i)
            hold on
            p1=plot(R.GRFs_separate(:,i),'-.','Color',Cs,'DisplayName','calcaneus');
            p2=plot(R.GRFs_separate(:,i+3),'--','Color',Cs,'DisplayName','forefoot');
            p3=plot(R.GRFs_separate(:,i+6),':','Color',Cs,'DisplayName','toes');
            title(R.colheaders.GRF{i});
            xlabel('% stride');
            ylabel('% body weight')
            
        end
%         if boolFirst
%             lh=legend([l,p1,p2,p3],'location','best');
%             lh.Interpreter = 'none';
%         end
        

        % get GRFs
        GRF1 = R.GRFs_separate(:,2);
        GRF2 = R.GRFs_separate(:,5);
        GRF3 = R.GRFs_separate(:,8);
        % total
        GRFt = GRF1 + GRF2 + GRF3;
        % relative
        grf1 = GRF1./GRFt*100;
        grf2 = GRF2./GRFt*100;
        grf3 = GRF3./GRFt*100;
        % set to 0 when total is too small
        grf1(abs(GRFt)<1) = 0;
        grf2(abs(GRFt)<1) = 0;
        grf3(abs(GRFt)<1) = 0;

        subplot(3,4,4)
        hold on
        p1=plot(x,grf1,'-.','Color',Cs,'DisplayName','calcaneus');
        p2=plot(x,grf2,'--','Color',Cs,'DisplayName','forefoot');
        p3=plot(x,grf3,':','Color',Cs,'DisplayName','toes');
        title(R.colheaders.GRF{2});
        xlabel('% stride');
        ylabel('% total GRF')
        
        lh=legend([p1,p2,p3],'location','best');
        lh.Interpreter = 'none';
        
        subplot(3,4,8)
        hold on
        plot(x,R.Muscle.vM(:,iSol),'Color',Cs)
        xlabel('% stride');
        ylabel('velocity')
        title('Soleus fibre velocity')
        grid on
            
        if isfield(R,'AnkleInGround')
            ictt = find(R.AnkleInGround.leverArmGRF.r~=0);
            relPos = (R.COPR - R.AnkleInGround.position.r)*1e3;
            mks = 3;
            
            subplot(3,4,5)
            hold on
            plot(x(ictt),relPos(ictt,1),'.','MarkerSize',mks,'Color',Cs)
            xlim([x(1),x(end)])
            xlabel('% stride');
            ylabel('distance (mm)')
            title('fore-aft distance COP wrt ankle')
            
            subplot(3,4,6)
            hold on
            plot(x(ictt),relPos(ictt,2),'.','MarkerSize',mks,'Color',Cs)
            xlim([x(1),x(end)])
            xlabel('% stride');
            ylabel('distance (mm)')
            title('vertical distance COP wrt ankle')
            
            subplot(3,4,7)
            hold on
            plot(x(ictt),relPos(ictt,3),'.','MarkerSize',mks,'Color',Cs)
            xlim([x(1),x(end)])
            xlabel('% stride');
            ylabel('distance (mm)')
            title('lateral distance COP wrt ankle')
            
            for ind=1:length(x)
                vec_a = R.COPR(ind,:) - R.AnkleInGround.position.r(ind,:); % from ankle origin to COP
                vec_b = R.GRFs(ind,1:3); 
                vec_ap = dot(vec_a,vec_b)/dot(vec_b,vec_b)*vec_b; % orthogonal projection of a onto b
                vec_an = vec_a - vec_ap; % component of a that is normal to b
                MA1(ind) = norm(vec_an);
                MA2(ind) = norm(vec_a);
                
                vec_a = vec_an;
                vec_b = R.AnkleInGround.axisOrientation.r(ind,:);
                vec_ap = dot(vec_a,vec_b)/dot(vec_b,vec_b)*vec_b; % orthogonal projection of a onto b
                vec_an = vec_a - vec_ap; % component of a that is normal to b
                MA(ind) = norm(vec_an);
                
                
            end
            
            subplot(3,4,9)
            hold on
            plot(x(ictt),MA(ictt)*1e3,'.','MarkerSize',mks,'Color',Cs)
            xlim([x(1),x(end)])
            xlabel('% stride');
            ylabel('lever arm (mm)')
            title('GRF lever arm wrt ankle')
            
            subplot(3,4,10)
            hold on
            plot(x(ictt),-R.AnkleInGround.leverArmSol.r(ictt)*1e3,'.','MarkerSize',mks,'Color',Cs)
            xlim([x(1),x(end)])
            xlabel('% stride');
            ylabel('lever arm (mm)')
            title('Soleus lever arm wrt ankle')
            
            subplot(3,4,11)
            hold on
            gearRatio = -MA(ictt)'./R.AnkleInGround.leverArmSol.r(ictt);
            plot(x(ictt),gearRatio,'.','MarkerSize',mks,'Color',Cs,'DisplayName',LegName);
            xlim([x(1),x(end)])
            xlabel('% stride');
            ylabel('Gear ratio (-)')
            title('GRF lever arm / sol lever')
            if boolFirst
                lh=legend('-DynamicLegend','location','east');
                lh.Interpreter = 'none';
            end
        end
    end

    %% Objective function
    axes('parent', tab8);
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
            lh=legend('-DynamicLegend','location','west');
            lh.Interpreter = 'none';
            lhPos = lh.Position;
            lhPos(1) = lhPos(1)+0.1;
            set(lh,'position',lhPos);
        end
    end
    
    
    %% Plot ankle muscle energetics
    axes('parent', tab9);
    
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
    
    if boolFirst && md 
        
        if~strcmp(subject,'Fal_s1')
            iSol_data = find(strcmp(Dat.(type).EMGheaders,'soleus_r'));
            iGas_data = find(strcmp(Dat.(type).EMGheaders,'gas_med_r'));
            iGas2_data = find(strcmp(Dat.(type).EMGheaders,'gas_lat_r'));
            iTib_data = find(strcmp(Dat.(type).EMGheaders,'tib_ant_r'));
            
            ankle_act(:,1) = Dat.(type).gc.lowEMG_mean(:,iSol_data);
            ankle_act(:,2) = Dat.(type).gc.lowEMG_mean(:,iGas_data);
            ankle_act(:,3) = Dat.(type).gc.lowEMG_mean(:,iGas2_data);
            ankle_act(:,4) = Dat.(type).gc.lowEMG_mean(:,iTib_data);
        
        else
            iSol_data = find(strcmp(Dat.(type).EMGheaders,'Soleus'));
            iGas_data = find(strcmp(Dat.(type).EMGheaders,'Gastrocnemius-medialis'));
            iGas2_data = find(strcmp(Dat.(type).EMGheaders,'Gastrocnemius-lateralis'));
            iTib_data = find(strcmp(Dat.(type).EMGheaders,'Tibialis-anterior'));
            
            ankle_act(:,1) = Dat.(type).gc.lowEMG_mean([51:end,1:50],iSol_data);
            ankle_act(:,2) = Dat.(type).gc.lowEMG_mean([51:end,1:50],iGas_data);
            ankle_act(:,3) = Dat.(type).gc.lowEMG_mean([51:end,1:50],iGas2_data);
            ankle_act(:,4) = Dat.(type).gc.lowEMG_mean([51:end,1:50],iTib_data);
            
        end
        
        
        
        ankle_a = [ankle_act(ceil(end/2):end,:); ankle_act(1:ceil(end/2)-1,:)];
        
%         sol_act_std = Dat.(type).gc.lowEMG_std(:,iSol_data);
%         gas_act_std = Dat.(type).gc.lowEMG_std(:,iGas_data);
%         gas2_act_std = Dat.(type).gc.lowEMG_std(:,iGas2_data);
%         tib_act_std = Dat.(type).gc.lowEMG_std(:,iTib_data);

        scale(1) = max(R.a(:,iSol))/max(ankle_act(:,1));
        scale(2) = max(R.a(:,iGas))/max(ankle_act(:,2));
        scale(3) = max(R.a(:,iGas2))/max(ankle_act(:,3));
        scale(4) = max(R.a(:,iTib))/max(ankle_act(:,4));
        
        scale = diag(scale);
        
        ankle_a_sc = ankle_a*scale;
        
        for i=1:4
            subplot(5,4,i)
            hold on
            plot(ankle_a_sc(:,i),'-k') 
        end

    end
    
    
    mVect = {'Soleus','Gas-lat','Gas-med','Tib-ant'};
    
    iM = [iSol iGas iGas2 iTib];
    
    for i=1:4
        subplot(5,4,i)
        hold on
        plot(R.a(:,iM(i)),'-','Color',Cs);  title(mVect{i});
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
    
    axes('parent', tab11);
    
    subplot(2,3,1); hold on;
    if boolFirst && md && ~strcmp(subject,'Fal_s1')
        plot(0,Dat.Passive.speed./Dat.Normal.Step.StrideFreq_mean,'o','Color','k','MarkerFaceColor','k','DisplayName','Data normal shoes');
        plot(0,Dat.Normal.speed./Dat.Normal.Step.StrideFreq_mean,'o','Color',darkgray,'MarkerFaceColor',darkgray,'DisplayName','Data normal shoes');
        plot(1,Dat.Active.speed./Dat.Active.Step.StrideFreq_mean,'o','Color',lightgray,'MarkerFaceColor',lightgray,'DisplayName','Data exo active');
    end
    plot(xParam,R.StrideLength,'o','Color',Cs,'MarkerFaceColor',Cs)
    ylabel('Stride length [m]');  xlabel(xParamLab);
    xlim([-0.5,1.5])
    
    subplot(2,3,2); hold on;
    if boolFirst && md && ~strcmp(subject,'Fal_s1')
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
    if boolFirst && md && ~strcmp(subject,'Fal_s1')
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
    if boolFirst && md && ~strcmp(subject,'Fal_s1')
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
    if boolFirst && md && ~strcmp(subject,'Fal_s1')
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
    if boolFirst && md && ~strcmp(subject,'Fal_s1')
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
    
    
    %% Windlass mechanism
    
    
    has_tmt = isfield(R.S,'tmt') && ~isempty(R.S.tmt) && R.S.tmt;
    has_tmt_unlocked =  isfield(R.S,'tmt_locked') && ~isempty(R.S.tmt_locked) && ~R.S.tmt_locked;
    has_WL = isfield(R.S,'Windlass') && ~isempty(R.S.Windlass) && R.S.Windlass ~= 0;
    
    has_mtj = isfield(R.S,'mtj') && ~isempty(R.S.mtj) && R.S.mtj;
    
    imtp = find(strcmp(R.colheaders.joints,'mtp_angle_r'));
    x = 1:(100-1)/(size(R.Qs,1)-1):100;
    
    if has_mtj
        imtj = find(strcmp(R.colheaders.joints,'mtj_angle_r'));
    else
        imtj = find(strcmp(R.colheaders.joints,'tmt_angle_r'));
    end
        
    axes('parent', tab12);
    
    subplot(3,5,13)
    hold on
    plot(x,R.TPass(:,imtp),'-','color',Cs,'linewidth',line_linewidth,'DisplayName','Passive');
    title('Passive mtp torque')
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    ylabel('Torque (Nm)','Fontsize',label_fontsize);
        
    subplot(3,5,2)
    hold on
    plot(x,R.Qs(:,imtp),'color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
    title('Mtp angle')
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    ylabel('Angle (�)','Fontsize',label_fontsize);
    
    subplot(3,5,7)
    hold on
    plot(x,R.Tid(:,imtp),'color',Cs,'linewidth',line_linewidth)
    title('Total mtp torque')
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    ylabel('Torque (Nm)','Fontsize',label_fontsize);
    
    subplot(3,5,5)
    hold on
    p1=plot(R.Qs(:,imtp),R.Tid(:,imtp),'-','color',Cs,'linewidth',line_linewidth,'DisplayName','total');
%     p2=plot(R.Qs(:,imtp),R.TPass(:,imtp),'--','color',Cs,'linewidth',line_linewidth,'DisplayName','Passive');
    xlabel('mtp angle (�)','Fontsize',label_fontsize);
    ylabel('mtp torque (Nm)','Fontsize',label_fontsize);
    title('mtp stiffness')
%     legend([p1,p2]);
    
%     subplot(3,5,15)
%     qk = zeros(length(x),1);
%     qs = [R.Qs(end-2:end,imtp); R.Qs(:,imtp); R.Qs(1:2,imtp)];
% %     ts = [R.TPass(end-2:end,imtp); R.TPass(:,imtp); R.TPass(1:2,imtp)];
%     ts = [R.Tid(end-2:end,imtp); R.Tid(:,imtp); R.Tid(1:2,imtp)];
%     for iii=1:length(x)
%         if R.Qs(iii,imtp)>pi/180
%             tmp = polyfit(qs(iii:iii+4),ts(iii:iii+4),1);
%             qk(iii) = -tmp(1);
%         end
%     end
%     hold on
%     plot(x,qk(:),'color',Cs,'linewidth',line_linewidth)
%     xlabel('mtp angle (�)','Fontsize',label_fontsize);
%     ylabel('mtp quasi-stiffness (Nm/rad)','Fontsize',label_fontsize);
%     title('mtp quasi-stiffness')

    
%     subplot(3,5,10)
%     hold on
%     
%     xlabel('mtp angle (�)','Fontsize',label_fontsize);
%     ylabel('mtp torque (Nm)','Fontsize',label_fontsize);
%     title('passive mtp stiffness')
    
    if isfield(R,'GRFs_separate') && ~isempty(R.GRFs_separate)
        subplot(3,5,14)
        hold on
        p1=plot(R.GRFs_separate(:,2),'-.','Color',Cs,'DisplayName','calcaneus');
        p2=plot(R.GRFs_separate(:,2+3),'--','Color',Cs,'DisplayName','forefoot');
        p3=plot(R.GRFs_separate(:,2+6),':','Color',Cs,'DisplayName','toes');
        title(R.colheaders.GRF{i});
        xlabel('% stride');
        ylabel('% body weight')
        title('Vertical GRF')
        lh=legend([p1,p2,p3],'location','best');
        lh.Interpreter = 'none';
    end
    
    if has_tmt && has_tmt_unlocked && has_WL || has_mtj
        
        q_mtj = R.Qs(:,imtj);
        qdot_mtj = R.Qdots(:,imtj)*pi/180;
        q_mtp = R.Qs(:,imtp);

        if isfield(R,'windlass') && ~isempty(R.windlass)
            M_PF = R.windlass.M_PF;
            M_li = R.windlass.M_li;
            M_mtp = R.windlass.M_mtp;
            F_PF = R.windlass.F_PF;
            l_PF = R.windlass.l_PF;
            l_fa = R.windlass.l_fa;
            h_fa = R.windlass.h_fa;
            MA_PF = R.windlass.MA_PF;
            L0 = R.windlass.L0;
            H0 = R.windlass.H0;
            
        else
            % Some earlier versions did not analyse the windlass mechanism
            % results in f_LoadSim.
            kTMT_li = 1.5/(pi/180)/5;
            kTMT_PF = R.S.kTMT;
            dTMT = R.S.dTMT;
            cWL = R.S.cWL;

            M_PF = zeros(length(x),1);
            M_li = zeros(length(x),1);
            M_mtp = zeros(length(x),1);
            F_PF = zeros(length(x),1);
            l_PF = zeros(length(x),1);
            l_fa = zeros(length(x),1);
            h_fa = zeros(length(x),1);
            MA_PF = zeros(length(x),1);

            if has_mtj
                AddCasadiPaths();
                import casadi.*
                f_PF_stiffness = f_getPlantarFasciaStiffnessModelCasADiFunction(R.S.PF_stiffness);
            end

            for i=1:length(R.Qs)
                if has_mtj
                    [~, M_PFi,F_PFi,~,~,~,li,~,L0,hi,~,H0,~] = ...
                        getPassiveMtjMomentWindlass_v2(q_mtj(i)*pi/180,qdot_mtj(i),...
                        q_mtp(i)*pi/180,f_PF_stiffness);

                else
                    [~, M_PFi,F_PFi,~,~,li,~,L0,hi,~,H0,~] = ...
                        getPassiveTmtjMomentWindlass(q_mtj(i)*pi/180,qdot_mtj(i),...
                        q_mtp(i)*pi/180,kTMT_li,kTMT_PF,dTMT,R.S.subject,cWL);
                    M_PFi = -M_PFi;
                end

                M_PF(i) = M_PFi;
                l_fa(i) = li;
                h_fa(i) = hi;
                F_PF(i) = F_PFi;
            end
        end

        
    
        subplot(3,5,10)
        hold on
        plot(R.Qs(:,imtp),l_PF,'-','color',Cs,'linewidth',line_linewidth)
        xlabel('mtp angle (�)','Fontsize',label_fontsize);
        ylabel('PF length','Fontsize',label_fontsize);
        title('mtp stiffness')
        
        subplot(3,5,1)
        hold on
        plot(x,q_mtj,'color',Cs,'linewidth',line_linewidth)
        title('Midtarsal angle')
        xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        ylabel('Angle (�)','Fontsize',label_fontsize);

        subplot(3,5,3)
        hold on
        plot(x,l_fa*1000,'color',Cs,'linewidth',line_linewidth)
        title('Foot arch length')
        xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        ylabel('Length (mm)','Fontsize',label_fontsize);
        
        subplot(3,5,4)
        hold on
        plot(x,h_fa*1000,'color',Cs,'linewidth',line_linewidth)
        title('Foot arch height')
        xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        ylabel('Height (mm)','Fontsize',label_fontsize);
        
        subplot(3,5,6)
        hold on
        p1=plot(x,R.Tid(:,imtj),'color',Cs,'linewidth',line_linewidth,'DisplayName','Total');
%         p2=plot(x,M_PF,':','color',Cs,'linewidth',line_linewidth,'DisplayName','Plantar fascia');
%         legend([p1,p2],'location','best')
        title('Midtarsal torque')
        xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        ylabel('Torque (Nm)','Fontsize',label_fontsize);
        
        subplot(3,5,11)
        hold on
        p1=plot(x,M_PF,'-','color',Cs,'linewidth',line_linewidth,'DisplayName','Plantar fascia');
        p2=plot(x,M_li,'--','color',Cs,'linewidth',line_linewidth,'DisplayName','Passive (non-PF)');
        legend([p1,p2],'location','best')
        title('Midtarsal torque')
        xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        ylabel('Torque (Nm)','Fontsize',label_fontsize);

        
        subplot(3,5,12)
        hold on
        p1=plot(x,R.Tid(:,imtp)-R.TPass(:,imtp),':','color',Cs,'linewidth',line_linewidth,'DisplayName','Active');
        if isfield(R.S,'WL_T_mtp') && R.S.WL_T_mtp
            p2=plot(x,M_mtp,'-','color',Cs,'linewidth',line_linewidth,'DisplayName','Plantar fascia');
            p3=plot(x,R.TPass(:,imtp)-M_mtp,'--','color',Cs,'linewidth',line_linewidth,'DisplayName','Passive (non-PF)');
            legend([p1,p2,p3],'location','best')
        else
            p3=plot(x,R.TPass(:,imtp),'--','color',Cs,'linewidth',line_linewidth,'DisplayName','Passive');
            legend([p1,p3],'location','best')
        end
        title('Mtp torque')
        xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        ylabel('Torque (Nm)','Fontsize',label_fontsize);
        
        subplot(3,5,8)
        hold on
        if isfield(R.S,'PF_slack_length')
            ls = R.S.PF_slack_length;
            dl = l_PF-ls;
            PF_strain = (l_PF./ls-1)*100;
            plot(x,PF_strain,'color',Cs,'linewidth',line_linewidth)
        else
            dl = l_PF-0.17;
        end
        title('Plantar fascia strain')
        xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        ylabel('Nominal strain (%)','Fontsize',label_fontsize);
        
        subplot(3,5,9)
        hold on
%         plot(x,F_PF/(R.body_mass*9.81)*100,'color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
        plot(x,F_PF,'color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
        title('Plantar fascia force')
        xlabel('Gait cycle (%)','Fontsize',label_fontsize);
%         ylabel('Force/BW (%)','Fontsize',label_fontsize);
        ylabel('Force (N)','Fontsize',label_fontsize);
        
    
    
        %% Windlass (exta tab)
        
        W_PF = zeros(length(x),1);
        W_li = zeros(length(x),1);
        for ii=2:length(x)
           W_PF(ii) = trapz(dl(1:ii),-F_PF(1:ii)); % sign convention!
           W_li(ii) = trapz(q_mtj(1:ii)*pi/180,M_li(1:ii));
        end
        W_PF(2:end) = (W_PF(2:end)-W_PF(2))/R.body_mass;
        W_li(2:end) = (W_li(2:end)-W_li(2))/R.body_mass;
        
        dt = R.tf_step*2/length(x);
        v_PF = zeros(length(x),1);
        for ii=1:length(x)
            i2h = mod(ii+1,length(x))+1; % loop around at the end
            i1h = mod(ii+0,length(x))+1;
            i0h = mod(ii-1,length(x))+1;
            v_PF(ii) = -(l_PF(i2h) - 4*l_PF(i1h) + 3*l_PF(i0h))/(2*dt);
        end
        P_PF = -v_PF.*F_PF/R.body_mass;
        P_li = qdot_mtj.*M_li/R.body_mass;
        
        
    else
        axes('parent', tab12);
        
        subplot(3,5,11)
        hold on
        p1=plot(x,R.Tid(:,imtp)-R.TPass(:,imtp),':','color',Cs,'linewidth',line_linewidth,'DisplayName','Active');
        p3=plot(x,R.TPass(:,imtp),'--','color',Cs,'linewidth',line_linewidth,'DisplayName','Passive');
        legend([p1,p3],'location','best')
        title('Mtp torque')
        xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        ylabel('Torque (Nm)','Fontsize',label_fontsize);

        subplot(3,5,2)
        lh = legend('location','best');
        lh.Interpreter = 'none';
    
    end
    
    W_mtj = zeros(length(x),1);
    if has_tmt && has_tmt_unlocked || has_mtj

        for ii=2:length(x)
           W_mtj(ii) = trapz(R.Qs(1:ii,imtj)*pi/180,R.Tid(1:ii,imtj));
        end
        W_mtj(2:end) = (W_mtj(2:end)-W_mtj(2))/R.body_mass;        
        P_mtj = R.Qdots(:,imtj)*pi/180.*R.Tid(:,imtj)/R.body_mass;
    else
        P_mtj = W_mtj;
    end
    
    iankle = strcmp(R.colheaders.joints,'ankle_angle_r');
    isubt = strcmp(R.colheaders.joints,'subtalar_angle_r');

    W_ankle = zeros(length(x),1);
    W_subt = zeros(length(x),1);
    W_mtp = zeros(length(x),1);
    for ii=2:length(x)
       W_ankle(ii) = trapz(R.Qs(1:ii,iankle)*pi/180,R.Tid(1:ii,iankle));
       W_subt(ii) = trapz(R.Qs(1:ii,isubt)*pi/180,R.Tid(1:ii,isubt));
       W_mtp(ii) = trapz(R.Qs(1:ii,imtp)*pi/180,R.Tid(1:ii,imtp));
    end
    W_ankle(2:end) = (W_ankle(2:end)-W_ankle(2))/R.body_mass;
    W_subt(2:end) = (W_subt(2:end)-W_subt(2))/R.body_mass;
    W_mtp(2:end) = (W_mtp(2:end)-W_mtp(2))/R.body_mass;
    W_tot = W_mtj + W_ankle + W_subt + W_mtp;

    P_ankle = R.Qdots(:,iankle)*pi/180.*R.Tid(:,iankle)/R.body_mass;
    P_subt = R.Qdots(:,isubt)*pi/180.*R.Tid(:,isubt)/R.body_mass;
    P_ankle_pass = R.Qdots(:,iankle)*pi/180.*R.TPass(:,iankle)/R.body_mass;
    P_subt_pass = R.Qdots(:,isubt)*pi/180.*R.TPass(:,isubt)/R.body_mass;
    P_mtp = R.Qdots(:,imtp)*pi/180.*R.Tid(:,imtp)/R.body_mass;
    P_tot = P_mtj + P_ankle + P_subt + P_mtp;

    if isfield(R,'vT')
       P_T_Sol = -R.FT(:,iSol).*R.vT(:,iSol)/R.body_mass;
       P_T_Gas = -R.FT(:,iGas).*R.vT(:,iGas)/R.body_mass;
       P_T_Gas2 = -R.FT(:,iGas2).*R.vT(:,iGas2)/R.body_mass;
    end
    
    axes('parent', tab14);

    % Work
    if has_tmt && has_tmt_unlocked || has_mtj
        subplot(3,4,1)
        hold on
        grid on
        plot(x,W_mtj,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
        ylabel('Work (J/kg)')
        xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        title('Midtarsal')

        subplot(4,4,4)
        hold on
        grid on
        plot(R.S.kMT_li,W_tot(end),'o','Color',Cs,'MarkerFaceColor',Cs)
%         xlabel('mtj stiffness (Nm/rad)')
        ylabel('Work (J/kg)')
        title('Mechanical work foot (/gait cycle)')

        dist_trav = R.Qs(end,strcmp(R.colheaders.joints,'pelvis_tx')) - ...
            R.Qs(1,strcmp(R.colheaders.joints,'pelvis_tx'));
        
        subplot(4,4,8)
        hold on
        grid on
        plot(R.S.kMT_li,R.COT*dist_trav,'o','Color',Cs,'MarkerFaceColor',Cs)
%         xlabel('mtj stiffness (Nm/rad)')
        ylabel('Metab. Energy (J/kg)')
        title('Metabolic energy (/gait cycle)')

        subplot(4,4,12)
        hold on
        grid on
        plot(R.S.kMT_li,W_tot(end)/dist_trav,'o','Color',Cs,'MarkerFaceColor',Cs)
%         xlabel('mtj stiffness (Nm/rad)')
        ylabel('Work (J/(kg m))')
        title({'Mechanical work foot (/m)'})
        
        subplot(4,4,16)
        hold on
        grid on
        plot(R.S.kMT_li,R.COT,'o','Color',Cs,'MarkerFaceColor',Cs)
        xlabel('mtj stiffness (Nm/rad)')
        ylabel('COT (J/(kg m)')
        title('Cost of transport')
    
        if has_WL || has_mtj
            subplot(3,4,5)
            hold on
            grid on
            plot(x,W_PF,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
            ylabel('Work (J/kg)')
            xlabel('Gait cycle (%)','Fontsize',label_fontsize);
            title('Plantar fascia')

            subplot(3,4,9)
            hold on
            grid on
            plot(x,W_li,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
            ylabel('Work (J/kg)')
            xlabel('Gait cycle (%)','Fontsize',label_fontsize);
            title('Other mtj')
        end
    end
    
    subplot(3,4,2)
    hold on
    grid on
    plot(x,W_ankle,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
    ylabel('Work (J/kg)')
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    title('Ankle')

    subplot(3,4,6)
    hold on
    grid on
    plot(x,W_subt,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
    ylabel('Work (J/kg)')
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    title('Subtalar')

    subplot(3,4,10)
    hold on
    grid on
    plot(x,W_mtp,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
    ylabel('Work (J/kg)')
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    title('Mtp')
    
    subplot(3,4,11)
    hold on
    grid on
    plot(x,W_mtp+W_mtj,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
    ylabel('Work (J/kg)')
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    title('Mtp + mtj')

    if boolFirst
        lh=legend('-DynamicLegend','location','west');
        lh.Interpreter = 'none';
        lhPos = lh.Position;
        lhPos(1) = lhPos(1)+0.1;
        set(lh,'position',lhPos);
    end

    subplot(3,4,3)
    hold on
    grid on
    plot(x,W_tot,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
    ylabel('Work (J/kg)')
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    title('Total foot')

    if isfield(R,'GRFs_separate') && ~isempty(R.GRFs_separate)

        subplot(3,4,7)
        hold on
        grid on
        p1=plot(R.GRFs_separate(:,2),'-.','Color',Cs,'DisplayName','calcaneus');
        p2=plot(R.GRFs_separate(:,2+3),'--','Color',Cs,'DisplayName','forefoot');
        p3=plot(R.GRFs_separate(:,2+6),':','Color',Cs,'DisplayName','toes');
        title(R.colheaders.GRF{i});
        xlabel('% stride');
        ylabel('% body weight')
        title('Vertical GRF')
        lh=legend([p1,p2,p3],'location','best');
        lh.Interpreter = 'none';
    end

    axes('parent', tab13);
    % Power
    if has_tmt && has_tmt_unlocked || has_mtj
        subplot(3,4,1)
        hold on
        grid on
        plot(x,P_mtj,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
        ylabel('Power (W/kg)')
        xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        title('Midtarsal')

        if has_WL || has_mtj
            subplot(3,4,5)
            hold on
            grid on
            plot(x,P_PF,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
            ylabel('Power (W/kg)')
            xlabel('Gait cycle (%)','Fontsize',label_fontsize);
            title('Plantar fascia')

            subplot(3,4,9)
            hold on
            grid on
            plot(x,P_li,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
            ylabel('Power (W/kg)')
            xlabel('Gait cycle (%)','Fontsize',label_fontsize);
            title('Other mtj')
        end
    end
    
    subplot(3,4,2)
    hold on
    grid on
    plot(x,P_ankle,'Color',Cs,'linewidth',line_linewidth,'DisplayName','Ankle');
    ylabel('Power (W/kg)')
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    title('Ankle')
    
    subplot(3,4,6)
    hold on
    grid on
    plot(x,P_subt,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
    ylabel('Power (W/kg)')
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    title('Subtalar')

    subplot(3,4,10)
    hold on
    grid on
    plot(x,P_mtp,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
    ylabel('Power (W/kg)')
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    title('Mtp')

    subplot(3,4,11)
    hold on
    grid on
    plot(x,P_mtp+P_mtj,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
    ylabel('Power (W/kg)')
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    title('Mtp + mtj')
    
    if boolFirst
        lh=legend('-DynamicLegend','location','west');
        lh.Interpreter = 'none';
        lhPos = lh.Position;
        lhPos(1) = lhPos(1)+0.1;
        set(lh,'position',lhPos);
    end

    subplot(3,4,3)
    hold on
    grid on
    plot(x,P_tot,'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
    ylabel('Power (W/kg)')
    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    title('Total foot')

    if isfield(R,'vT')
        subplot(3,4,4)
        hold on
        grid on
        plot(x,P_T_Sol,'-','Color',Cs,'linewidth',line_linewidth)
        ylabel('Power (W/kg)')
        xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        title('Tendon Soleus')
    
        subplot(3,4,8)
        hold on
        grid on
%         p1=plot(x,P_T_Gas,'-','Color',Cs,'linewidth',line_linewidth,'DisplayName','Gas lat');
%         p2=plot(x,P_T_Gas2,':','Color',Cs,'linewidth',line_linewidth,'DisplayName','Gas med');
        plot(x,P_T_Gas+P_T_Gas2,'-','Color',Cs,'linewidth',line_linewidth)
        ylabel('Power (W/kg)')
        xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        title('Tendon Gastrocnemius')
%         legend([p1,p2]);
    end
    
    if isfield(R,'GRFs_separate') && ~isempty(R.GRFs_separate)
        subplot(3,4,7)
        hold on
        grid on
        p1=plot(R.GRFs_separate(:,2),'-.','Color',Cs,'DisplayName','calcaneus');
        p2=plot(R.GRFs_separate(:,2+3),'--','Color',Cs,'DisplayName','forefoot');
        p3=plot(R.GRFs_separate(:,2+6),':','Color',Cs,'DisplayName','toes');
        title(R.colheaders.GRF{i});
        xlabel('% stride');
        ylabel('% body weight')
        title('Vertical GRF')
        lh=legend([p1,p2,p3],'location','best');
        lh.Interpreter = 'none';
    end
        
        
    
   
    %% Exo assistance

    if isfield(R,'w_RotVel_exo') && ~isempty(R.w_RotVel_exo)
        axes('parent', tab15);
        
        T_mean = mean(R.T_exo(:,2));
        subplot(2,2,1)
        hold on
        plot(R.T_exo(:,2),'-','Color',Cs)
        line(get(gca, 'xlim'),[1,1]*T_mean,'color',Cs,'LineStyle','--')
        if R.S.T_max_ankle_exo > 0
            line(get(gca, 'xlim'),-[1,1]*R.S.T_max_ankle_exo,'color',Cs,'LineStyle',':')
            line(get(gca, 'xlim'),-[1,1]*R.S.T_min_ankle_exo,'color',Cs,'LineStyle',':')
        end
        ylabel('Exo Moment [Nm]')
        xlabel('% stride')
        title('Exo Moment')
        
        subplot(2,2,2)
        hold on
        plot(R.w_RotVel_exo(:,2)*180/pi,'-','Color',Cs)
        ylabel('Exo Velocity [�/s]')
        xlabel('% stride')
        title('Exo rotational velocity')
        
        P_mean = mean(R.P_exo(:,2));
        subplot(2,2,3)
        hold on
        plot(R.P_exo(:,2),'-','Color',Cs)
        line(get(gca, 'xlim'),[1,1]*P_mean,'color',Cs,'LineStyle','--')
        if R.S.P_max_ankle_exo > 0
            line(get(gca, 'xlim'),[1,1]*R.S.P_max_ankle_exo,'color',Cs,'LineStyle',':')
        end
        ylabel('Exo Power [W]')
        xlabel('% stride')
        title('Exo Power')
        
        subplot(2,2,4)
        title('Energy')
        hold on
        yyaxis left
        plot(1,sum(R.W_exo_rel),'o','Color',Cs)
        ylabel('Exo work (J/(kg m)')
        
        yyaxis right
        plot(2,R.COT,'o','Color',Cs)
        ylabel('COT (J/(kg m)')
        xlim([0,3])
    
    end
    
    
    %% Plot mtp muscle energetics
    axes('parent', tab10);
    
    ifd = find(strcmp(R.colheaders.muscles,'flex_dig_r'));
    ifh = find(strcmp(R.colheaders.muscles,'flex_hal_r'));
    ied = find(strcmp(R.colheaders.muscles,'ext_dig_r'));
    ieh = find(strcmp(R.colheaders.muscles,'ext_hal_r'));
    
    mVect = {'Flex-dig','Flex-hal','Ext-dig','Ext-hal'};
    
    iM = [ifd ifh ied ieh];
    
    for i=1:4
        subplot(5,4,i)
        hold on
        plot(R.a(:,iM(i)),'-','Color',Cs);  title(mVect{i});
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
    
    
else
    warning(['File not found: ' ResultsFile]);
end

end

