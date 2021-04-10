function [] = PlotResults_3DSim_Report(ResultsFile,LegNames,RefData,mtj,makeplot,figNamePrefix)

if strcmp(RefData,'none')
    md = 0;
else
    md = 1;
    if strcmp(RefData,'act')
        type = 'Active';
        refName = 'Powered exoskeleton';
    elseif strcmp(RefData,'pas')
        type = 'Passive';
        refName = 'Unpowered exoskeleton';
    elseif strcmp(RefData,'norm')
        type = 'Normal';
        refName = 'Shoes';
    elseif strcmp(RefData,'Fal_s1')
        type = 'Normal';
        refName = 'Barefoot';
    end
end

nr = length(ResultsFile);

if numel(LegNames) == nr
    LN = 1;
else
    LN = 0;
end

if strcmp(RefData,'Fal_s1')
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
        'Midtarsal L','Midtarsal R','Tarsometatarsal L','Tarsometatarsal R',...
        'Metatarsophalangeal L','Metatarsophalangeal R',...
        'Lumbar extension','Lumbar bending','Lumbar rotation',...
        'Arm flexion L','Arm adduction L','Arm rotation L',...
        'Arm flexion R','Arm adduction R','Arm rotation R',...
        'Elbow flexion L','Elbow flexion R'};
        
% prep figures
scs = get(0,'ScreenSize');
fwide = [scs(3)/2, scs(4)/2-100];
fhigh = [scs(3)/2, scs(4)-120];
flong = [scs(3)/2, scs(4)/4];

fpos = [1,scs(4)/2+20;
        1,40;
        -scs(3)/2,40;
        -scs(3),40;];


for inr=1:nr
    load(ResultsFile{inr},'R');
    if LN
        LegName = LegNames{inr};
    else
        LegName = R.S.savename;
    end
    
    has_no_tmt = ~isfield(R.S,'tmt') || isempty(R.S.tmt) || R.S.tmt == 0;
    has_no_mtj = ~isfield(R.S,'mtj') || isempty(R.S.mtj) || R.S.mtj == 0;
    
    if inr==1 && md
        pc_name = getenv('COMPUTERNAME');
        if strcmp(pc_name,'MSI')
            if strcmp(RefData,'Fal_s1')
                load('D:\school\WTK\thesis\model\3dpredictsim\Data\Fal_s1.mat','Dat');
            else
                % load data Pog_s1 from struct saved during ...\Analyze_ExoData\Batch\BatchScript_LatexReport.m
                load('D:\school\WTK\thesis\model\3dpredictsim\Data\Pog_s1.mat','Dat');
            end
            Qref = Dat.(type).gc;
        else
            md = 0;
        end
        
    end
    
    idx_js = [4,5,7,8,9,10,12]; % colheaders mocap
    idx_sp = [1,2,3,5,6,7,8]; % subplots to use
    
    if has_no_tmt && has_no_mtj
        idx_Qs = [10,11,14,16,18,20]; % Qs in result file
    else
        idx_Qs = [10,11,14,16,18,20,22];
    end
    
    if mtj==1
        idx_title = [10,11,14,16,18,20,24]; % names for title
    elseif mtj==0
        idx_title = [10,11,14,16,18,22,24];
    else
        idx_title = [10,11,14,16,18,24];
        idx_js = [4,5,7,8,9,12]; % colheaders mocap
    idx_sp = [1,2,3,5,6,7]; % subplots to use
    end
    
    
    %% kinematics
    if makeplot.kinematics
    
        if inr==1
            h1 = figure('Position',[fpos(1,:),fwide]);
        else
            figure(h1);
        end
        j = 0;
        label_fontsize  = 12;
        line_linewidth  = 0.5;
        NumTicks = 6;
        for i = 1:length(idx_title)
            subplot(2,4,idx_sp(i))
            x = 1:(100-1)/(size(R.Qs,1)-1):100;
            % Experimental data
            if  inr == 1 && md
                idx_jref = strcmp(Qref.colheaders,joints_ref{idx_js(i)});
                if sum(idx_jref) == 1
                    meanPlusSTD = (Qref.Qall_mean(:,idx_jref) + 2*Qref.Qall_std(:,idx_jref)).*180/pi;
                    meanMinusSTD = (Qref.Qall_mean(:,idx_jref) - 2*Qref.Qall_std(:,idx_jref)).*180/pi;

                    stepQ = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
                    intervalQ = 1:stepQ:size(R.Qs,1);
                    sampleQ = 1:size(R.Qs,1);
                    meanPlusSTD = interp1(intervalQ,meanPlusSTD,sampleQ);
                    meanMinusSTD = interp1(intervalQ,meanMinusSTD,sampleQ);

                    hold on
                    fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k','DisplayName',['MoCap ' refName]);
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
                plot(x,R.Qs(:,idx_Qs(j)),'linewidth',line_linewidth,'DisplayName',LegName);

            end

            % Plot settings
            if inr == nr
                set(gca,'Fontsize',label_fontsize);
                title(joints_tit{idx_title(i)},'Fontsize',label_fontsize);
                % Y-axis
                if i == 1 || i == 4
                    ylabel('Angle (°)','Fontsize',label_fontsize);
                end
                % X-axis
                L = get(gca,'XLim');
                set(gca,'XTick',linspace(L(1),L(2),NumTicks))
                if i > 3
                    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
                end
            end
            if i == 3
                lh=legend('-DynamicLegend','location','west');
                lh.Interpreter = 'none';
                lhPos = lh.Position;
                lhPos(1) = lhPos(1)+0.2;
    %             lhPos(2) = lhPos(2)-0.2;
                set(lh,'position',lhPos);
            end
        end
    end
    
    
    %% kinetics
    if makeplot.kinetics
        
        if inr==1
            h2 = figure('Position',[fpos(2,:),fwide]);
        else
            figure(h2);
        end
        j = 0;
        label_fontsize  = 12;
        line_linewidth  = 0.5;
        for i = 1:length(idx_title)
            subplot(2,4,idx_sp(i))
            x = 1:(100-1)/(size(R.Qs,1)-1):100;
            % Experimental data
            if  inr == 1 && md
                idx_jref = strcmp(Qref.colheaders,joints_ref{idx_js(i)});
                if sum(idx_jref) == 1
                    meanPlusSTD = Qref.Tall_mean(:,idx_jref) + 2*Qref.Tall_std(:,idx_jref);
                    meanMinusSTD = Qref.Tall_mean(:,idx_jref) - 2*Qref.Tall_std(:,idx_jref);
                    stepID = (size(R.Qs,1)-1)/(size(meanPlusSTD,1)-1);
                    intervalID = 1:stepID:size(R.Qs,1);
                    sampleID = 1:size(R.Qs,1);
                    meanPlusSTD = interp1(intervalID,meanPlusSTD,sampleID);
                    meanMinusSTD = interp1(intervalID,meanMinusSTD,sampleID);
                    hold on
                    fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k','DisplayName',['MoCap ' refName]);
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
                plot(x,R.Tid(:,idx_Qs(j)),'linewidth',line_linewidth,'DisplayName',LegName);

            end
            % Plot settings
            if inr == nr
                % Plot settings
                set(gca,'Fontsize',label_fontsize);
                title(joints_tit{idx_title(i)},'Fontsize',label_fontsize);
                % Y-axis
                if i == 1 || i == 4
                    ylabel('Torque (Nm)','Fontsize',label_fontsize);
                end
                % X-axis
                L = get(gca,'XLim');
                set(gca,'XTick',linspace(L(1),L(2),NumTicks))
                if i > 3
                    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
                end

            end
            if i == 3
                lh=legend('-DynamicLegend','location','west');
                lh.Interpreter = 'none';
                lhPos = lh.Position;
                lhPos(1) = lhPos(1)+0.2;
    %             lhPos(2) = lhPos(2)-0.2;
                set(lh,'position',lhPos);
            end
        end
    end
    
    
    %% soleus
    if makeplot.soleus
        
        if inr==1
            h3 = figure('Position',[fpos(3,:),fhigh]);
        else
            figure(h3);
        end
        iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));

        subplot(5,1,1); hold on;

        if inr==1 && md && ~strcmp(RefData,'Fal_s1')
            iSol_data = find(strcmp(Dat.(type).EMGheaders,'soleus_r'));
            ankle_act(:,1) = Dat.(type).gc.lowEMG_mean(:,iSol_data);
            scale = max(R.a(:,iSol))/max(ankle_act(:,1));
            ankle_a_sc = ankle_a.*scale;
            plot(ankle_a_sc(:,1),'-k','DisplayName','EMG data') 
        end

        plot(R.a(:,iSol),'-') 
        title('Soleus')
        ylabel('activity')

        subplot(5,1,2)
        plot(R.MetabB.Etot(:,iSol),'-'); hold on;
        ylabel('Muscle metab power');

        subplot(5,1,3)
        plot(R.lMtilde(:,iSol),'-'); hold on;
        ylabel('Norm fiber length');

        subplot(5,1,4)
        plot(R.MetabB.Wdot(:,iSol),'-'); hold on;
        ylabel('Wdot');

        subplot(5,1,5)
        plot(R.FT(:,iSol),'-','DisplayName',LegName); hold on;
        ylabel('Norm muscle force');
        xlabel('Gait cycle (%)','Fontsize',label_fontsize);


        if inr==1
            lh=legend('-DynamicLegend','location','east');
            lh.Interpreter = 'none';
        end
    end
    
    %% GRF
    if makeplot.GRF
        
        if inr==1
            h4 = figure('Position',[fpos(4,:),flong]);
        else
            figure(h4);
        end
        if inr==1 && md && ~strcmp(RefData,'Fal_s1')
            for i=1:3
                subplot(1,3,i)
                hold on
                plot(Dat.(type).gc.GRF.Fmean(:,i)/(R.body_mass*9.81)*100,'-k');
            end

        end
        for i=1:3
            subplot(1,3,i)
            hold on
            l = plot(R.GRFs(:,i),'-');
            title(R.colheaders.GRF{i});
            xlabel('Gait cycle (%)','Fontsize',label_fontsize);
            if i==1
                ylabel({'Ground reaction force','(% body weight)'})
            end
        end
        l.DisplayName = LegName;

        if inr == 1
            lh=legend('location','northeast');
            lh.Interpreter = 'none';
        end
    end

%     if isfield(R,'GRFs_separate') && ~isempty(R.GRFs_separate)
%         for i=1:3
%             subplot(1,3,i)
%             hold on
%             p1=plot(R.GRFs_separate(:,i),'-.','DisplayName','calcaneus');
%             p2=plot(R.GRFs_separate(:,i+3),'--','DisplayName','forefoot');
%             p3=plot(R.GRFs_separate(:,i+6),':','DisplayName','toes');
%             title(R.colheaders.GRF{i});
%             xlabel('% stride');
%             ylabel('% body weight')
%             
%         end
%         if boolFirst
%             lh=legend([l,p1,p2,p3],'location','best');
%             lh.Interpreter = 'none';
%         end
%     end
    
    
    %% windlass
    
    
end


if figNamePrefix ~= 0
    
    set(h1,'PaperPositionMode','auto')
    print(h1,[figNamePrefix '_qs'],'-dpng','-r0')
    set(h2,'PaperPositionMode','auto')
    print(h2,[figNamePrefix '_Ts'],'-dpng','-r0')
    set(h3,'PaperPositionMode','auto')
    print(h3,[figNamePrefix '_sol'],'-dpng','-r0')
    set(h4,'PaperPositionMode','auto')
    print(h4,[figNamePrefix '_GRF'],'-dpng','-r0')
end

end


