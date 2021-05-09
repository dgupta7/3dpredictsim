function [] = PlotResults_3DSim_Report(ResultsFile,LegNames,RefData,mtj,makeplot,figNamePrefix)

set(0,'defaultTextInterpreter','tex');

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

label_fontsize  = 12;
line_linewidth  = 0.5;
NumTicks = 6;
CsV = hsv(nr);
   
        
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
    
    x = 1:(100-1)/(size(R.Qs,1)-1):100;
    
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
                plot(x,R.Qs(:,idx_Qs(j)),'linewidth',line_linewidth,'Color',CsV(inr,:),'DisplayName',LegName);

            end

            % Plot settings
            if inr == nr
                set(gca,'Fontsize',label_fontsize);
                ttl_tmp = joints_tit{idx_title(i)};
                title(ttl_tmp(1:end-2),'Fontsize',label_fontsize);
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
                plot(x,R.Tid(:,idx_Qs(j)),'linewidth',line_linewidth,'Color',CsV(inr,:),'DisplayName',LegName);

            end
            % Plot settings
            if inr == nr
                % Plot settings
                set(gca,'Fontsize',label_fontsize);
                ttl_tmp = joints_tit{idx_title(i)};
                title(ttl_tmp(1:end-2),'Fontsize',label_fontsize);
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
        iGas = find(strcmp(R.colheaders.muscles,'lat_gas_r'));
        iGas2 = find(strcmp(R.colheaders.muscles,'med_gas_r'));
        if isempty(iGas)
            iGas = find(strcmp(R.colheaders.muscles,'gaslat_r'));
        end
        if isempty(iGas2)
            iGas2 = find(strcmp(R.colheaders.muscles,'gasmed_r'));
        end
        imus = [iSol,iGas,iGas2];
        
        
%         figure
%         subplot(211)
%         plot(R.a(:,iSol))
%         subplot(212)
%         plot(R.e(:,iSol))
        
        if inr==1 && md 
            if strcmp(RefData,'Fal_s1')
                iSol_data = find(strcmp(Dat.(type).EMGheaders,'Soleus'));
                iGas_data = find(strcmp(Dat.(type).EMGheaders,'Gastrocnemius-medialis'));
                iGas2_data = find(strcmp(Dat.(type).EMGheaders,'Gastrocnemius-lateralis'));
                
                ankle_act(:,1) = Dat.(type).gc.lowEMG_mean([51:end,1:50],iSol_data);
                ankle_act(:,2) = Dat.(type).gc.lowEMG_mean([51:end,1:50],iGas_data);
                ankle_act(:,3) = Dat.(type).gc.lowEMG_mean([51:end,1:50],iGas2_data);
                
            else
                iSol_data = find(strcmp(Dat.(type).EMGheaders,'soleus_r'));
                iGas_data = find(strcmp(Dat.(type).EMGheaders,'gas_med_r'));
                iGas2_data = find(strcmp(Dat.(type).EMGheaders,'gas_lat_r'));
                ankle_act(:,1) = Dat.(type).gc.lowEMG_mean([51:end,1:50],iSol_data);
                ankle_act(:,2) = Dat.(type).gc.lowEMG_mean([51:end,1:50],iGas_data);
                ankle_act(:,3) = Dat.(type).gc.lowEMG_mean([51:end,1:50],iGas2_data);
            end
            ankle_a = [ankle_act(ceil(end/2):end,:); ankle_act(1:ceil(end/2)-1,:)];
            
            for imu=1:3
                subplot(6,3,imu); hold on;
                yyaxis right
                plot(ankle_a(:,imu),'-k','DisplayName','EMG data')
                a1 = gca;
                a1.YColor = [0,0,0];
                if imu==3
                    ylabel('EMG (mV)')
                end
                yyaxis left
                a1 = gca;
                a1.YColor = [0,0,0];
            end
            
        end
        
        subplot(6,3,1); hold on;
        plot(R.a(:,iSol),'-','Color',CsV(inr,:),'DisplayName',LegName);
        title('Soleus')
        ylabel('Activity (-)')
        
        subplot(6,3,2); hold on;
        plot(R.a(:,iSol),'-','Color',CsV(inr,:),'DisplayName',LegName);
        title('Gastrocnemius-medialis')
        
        subplot(6,3,3); hold on;
        plot(R.a(:,iSol),'-','Color',CsV(inr,:),'DisplayName',LegName);
        title('Gastrocnemius-lateralis')
        
        if inr==1
            lh=legend('-DynamicLegend','location','northeast');
            lh.Interpreter = 'none';
        end
        
        for imu=1:3
            
            subplot(6,3,3+imu)
            plot(R.FT(:,imus(imu)),'-','Color',CsV(inr,:),'DisplayName',LegName);
            hold on
            if imu==1
                ylabel('F_{normal,muscle} (N)');
            end

            subplot(6,3,6+imu)
            plot(R.MetabB.Etot(:,imus(imu)),'-','Color',CsV(inr,:)); hold on;
            hold on
            if imu==1
                ylabel('P_{metab, muscle} (W)');
            end

            subplot(6,3,9+imu)
            plot(R.MetabB.Wdot(:,imus(imu)),'-','Color',CsV(inr,:)); hold on;
            hold on
            if imu==1
                ylabel('P_{mech, MTU} (W)');
            end
            
            subplot(6,3,12+imu)
            plot(R.lMtilde(:,imus(imu)),'-','Color',CsV(inr,:)); hold on;
            hold on
            if imu==1
                ylabel('Norm fiber length (-)');
            end

            subplot(6,3,15+imu)
            plot(R.Muscle.vM(:,imus(imu)),'Color',CsV(inr,:),'DisplayName',LegName);
            hold on
            if imu==1
                ylabel('Fibre velocity (s^{-1})')
            end

            xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        end
    end
    
    lhPos = lh.Position;
    lhPos(2) = lhPos(2)+0.08;
    set(lh,'position',lhPos);
    
    %% GRF
    if makeplot.GRF
        GRF_title = {'fore-aft','vertical','lateral'};
        if inr==1
            h4 = figure('Position',[fpos(4,:),flong]);
        else
            figure(h4);
        end
        if inr==1 && md
            for i=1:3
                subplot(1,7,[2*i-1,2*i])
                hold on
                if strcmp(RefData,'Fal_s1')
                    p0=plot(Dat.(type).gc.GRF.Fmean(:,i),'-k','DisplayName','Total (experimental)');
                else
                    p0=plot(Dat.(type).gc.GRF.Fmean(:,i)/(R.body_mass*9.81)*100,'-k','Total (experimental)');
                end
            end

        end
        for i=1:3
            subplot(1,7,[2*i-1,2*i])
            hold on
            l = plot(R.GRFs(:,i),'-','Color',CsV(inr,:));
            title(GRF_title{i});
            xlabel('Gait cycle (%)','Fontsize',label_fontsize);
            if i==1
                ylabel({'Ground reaction force','(% body weight)'})
            end
        end
        l.DisplayName = ['Total (' LegName ')'];

        if isfield(R,'GRFs_separate') && ~isempty(R.GRFs_separate)
            for i=1:3
                subplot(1,7,[2*i-1,2*i])
                hold on
                p1=plot(R.GRFs_separate(:,i),'-.','Color',CsV(inr,:),'DisplayName','hindfoot');
                p2=plot(R.GRFs_separate(:,i+3),'--','Color',CsV(inr,:),'DisplayName','forefoot');
                p3=plot(R.GRFs_separate(:,i+6),':','Color',CsV(inr,:),'DisplayName','toes');
%                 title(R.colheaders.GRF{i});
%                 xlabel('% stride');
%                 ylabel('% body weight')

            end
            if inr == 1
                lh=legend([p0,l,p1,p2,p3],'location','northeast');
                lh.Interpreter = 'none';
                lhPos = lh.Position;
                lhPos(1) = lhPos(1)+0.2;
                set(lh,'position',lhPos);
            end
        elseif inr == 1
            lh=legend('location','northeast');
            lh.Interpreter = 'none';
            lhPos = lh.Position;
            lhPos(1) = lhPos(1)+0.2;
            set(lh,'position',lhPos);
        end
    end

    
    
    %% literature
    
    if makeplot.compareLiterature
        
        if inr==1
            pathmain = pwd;
            [pathRepo,~,~]  = fileparts(pathmain);
            folder = '\Figures';
            file = 'calcn_foreft_Caravaggi18.png';
            pathRefImg = fullfile(pathRepo,folder,file);
            img_q_mtj = imread(pathRefImg);
            file = 'mtp1_Caravaggi18.png';
            pathRefImg = fullfile(pathRepo,folder,file);
            img_q_mtp = imread(pathRefImg);
            file = 'PF_force_Caravaggi09.png';
            pathRefImg = fullfile(pathRepo,folder,file);
            img_F_PF = imread(pathRefImg);
            file = 'foot_power_Takahashi17.png';
            pathRefImg = fullfile(pathRepo,folder,file);
            img_P_foot = imread(pathRefImg);
            file = 'mtj_power_Takahashi17.png';
            pathRefImg = fullfile(pathRepo,folder,file);
            img_P_mtj = imread(pathRefImg);
            
            h5 = figure('Position',[fpos(1,:),fwide]);
            
            subplot(121)
            hold on
            axis tight
            hi1 = image([-12,102],flip([-72,27]+30),img_q_mtj);
            uistack(hi1,'bottom')
            xlabel('Gait cycle (%)')
            ylabel('Midtarsal angle (°)')
            title('Joint angles,')

            subplot(122)
            hold on
            axis tight
            hi2 = image([-12,102],flip([-36,77]),img_q_mtp);
            uistack(hi2,'bottom')
            xlabel('Gait cycle (%)')
            ylabel('Mtp angle (°)')
            title('Caravaggi (2018), experiments + ID')
            lh=legend('location','southwest');
            lh.Interpreter = 'none';
            
            h6 = figure('Position',[fpos(4,:),scs(3)/4*3, scs(4)/2]);
            
            subplot(121)
            hold on
            axis tight
            hi3 = image([-4,108],flip([-1.71,3.6]),img_P_foot);
            uistack(hi3,'bottom')
            ylabel('Power (W/kg)');
            xlabel('Stance phase (%)');
            title('Joint power,')
            lh=legend;
            lh.Position = [0.1539    0.6207    0.1132    0.2269];
            lh.Interpreter = 'none';
            title(lh,'Distal to')

            subplot(122)
            hold on
            axis tight
            hi4 = image([-6.5,109],flip([-1.64,1.69]),img_P_mtj);
            uistack(hi4,'bottom')
            ylabel('Power (W/kg)');
            xlabel('Stance phase (%)');
            title('Takahashi (2017), experiment-based')
            lh=legend;
            lh.Position=[0.6234    0.6373    0.1458    0.1287];
            lh.Interpreter = 'none';
            
            h7 = figure('Position',[fpos(2,:),scs(3)/2, scs(4)/3]);
            title('Plantar fascia force, Caravaggi (2009), FE-simulation')
            hold on
            axis tight
            hi5 = image([-13,112],flip([-0.16,1.75]),img_F_PF);
            uistack(hi5,'bottom')
            ylabel('Tension/BW (-)');
            xlabel('Stance phase (%)');
%             title('Plantar fascia force')
            lh=legend('location','northwest');
            lh.Interpreter = 'none';
            
        end
        
        istance = 1:1:round(R.Event.Stance)+2;
        xst = linspace(1,100,length(istance));
        imtp = find(strcmp(R.colheaders.joints,'mtp_angle_r'));
        imtj = find(strcmp(R.colheaders.joints,'mtj_angle_r'));
        iankle = strcmp(R.colheaders.joints,'ankle_angle_r');
        isubt = strcmp(R.colheaders.joints,'subtalar_angle_r');
        line_linewidth = 2;
        
        % Q
        figure(h5)
        if ~isempty(imtj)
            subplot(121)
            plot(x,R.Qs(:,imtj),'linewidth',line_linewidth,'Color',CsV(inr,:),'DisplayName',LegName);
       
        end
        
        subplot(122)
        plot(x,R.Qs(:,imtp),'linewidth',line_linewidth,'Color',CsV(inr,:),'DisplayName',LegName);
        
        % P
        figure(h6)
        
        P_ankle = R.Qdots(istance,iankle)*pi/180.*R.Tid(istance,iankle)/R.body_mass;
        P_subt = R.Qdots(istance,isubt)*pi/180.*R.Tid(istance,isubt)/R.body_mass;
        P_mtp = R.Qdots(istance,imtp)*pi/180.*R.Tid(istance,imtp)/R.body_mass;
        if ~isempty(imtj)
            P_mtj = R.Qdots(istance,imtj)*pi/180.*R.Tid(istance,imtj)/R.body_mass;
            P_hindfoot = P_mtj+P_mtp;
            P_tot = P_mtj + P_ankle + P_subt + P_mtp;
            
            subplot(121)
            plot(xst,P_mtp','-.','Color',[0.7, 0.1, 0.1],'linewidth',line_linewidth,'DisplayName',['forefoot (' LegName ')']);
            plot(xst,P_hindfoot,'-.','Color',[0, 0.4470, 0.7410],'linewidth',line_linewidth,'DisplayName',['hindfoot (' LegName ')']);
            plot(xst,P_tot,'-.k','linewidth',line_linewidth,'DisplayName',['shank (' LegName ')']);
        
        else
            P_hindfoot = P_mtp;
            P_tot = P_ankle + P_subt + P_mtp;
            
            subplot(121)
            plot(xst,P_mtp','-.','linewidth',line_linewidth,'DisplayName',['forefoot (' LegName ')']);
            plot(xst,P_hindfoot,'-.','linewidth',line_linewidth,'DisplayName',['hindfoot (' LegName ')']);
            plot(xst,P_tot,'-.','linewidth',line_linewidth,'DisplayName',['shank (' LegName ')']);
        
        end
        
        

        subplot(122)
        if ~isempty(imtj)
            plot(xst,P_mtp,'-.','Color',[0.7, 0.1, 0.1],'linewidth',line_linewidth,'DisplayName',['Distal to forefoot (' LegName ')']);
            plot(xst,P_mtj,'-','Color',[0.3010, 0.7450, 0.9330],'linewidth',line_linewidth,'DisplayName',['Midtarsal joint (' LegName ')']);
            plot(xst,P_hindfoot,'-.','Color',[0, 0.4470, 0.7410],'linewidth',line_linewidth,'DisplayName',['Distal to hindfoot (' LegName ')']);
        else
            plot(xst,P_hindfoot,'-.','linewidth',line_linewidth,'DisplayName',['Distal to hindfoot (' LegName ')']);
        end
        
        
        % F
        if isfield(R,'windlass') && ~isempty(R.windlass)
            figure(h7)
            plot(xst,R.windlass.F_PF(istance)/(R.body_mass*9.81),'linewidth',line_linewidth,'Color',CsV(inr,:),'DisplayName',LegName);
        end
            
        


    end
    
    
    %% centre of pressure
    if makeplot.COP
        if inr==1
            h8 = figure('Position',[fpos(2,:),fwide]);
            
        end
        
        if isfield(R,'AnkleInGround')
            ictt = find(R.AnkleInGround.leverArmGRF.r~=0);
            relPos = (R.COPR(ictt,:) - R.AnkleInGround.position.r(ictt,:))*1e3;
            
            figure(h8)
            subplot(121)
            p1=plot(0,0,'d','Color',CsV(inr,:),'DisplayName',['Ankle (' LegName ')']);
            hold on
            grid on
            axis equal
            xlabel('fore-after (mm)')
            ylabel('vertical (mm)')
            title('Centre of pressure during stance')
            plot(relPos(:,1),relPos(:,2),'o','Color',p1.Color,'DisplayName',['COP (' LegName ')']);
            
            subplot(122)
            p1=plot(0,0,'d','Color',CsV(inr,:),'DisplayName',['Ankle (' LegName ')']);
            hold on
            grid on
            axis equal
            ylabel('fore-after (mm)')
            xlabel('lateral (mm)')
            title('Centre of pressure during stance')
            plot(relPos(:,3),relPos(:,1),'o','Color',p1.Color,'DisplayName',['COP (' LegName ')']);
            legend('location','best')
            
            
%             figure
%             ictt1 = ictt(1:5:end);
%             grfx = [R.COPR(ictt1,1)*1e3, R.COPR(ictt1,1)*1e3 + R.GRFs(ictt1,1)];
%             grfy = [R.COPR(ictt1,2)*1e3, R.COPR(ictt1,2)*1e3 + R.GRFs(ictt1,2)];
%             plot(R.AnkleInGround.position.r(ictt1,1)*1e3,R.AnkleInGround.position.r(ictt1,2)*1e3,'d','DisplayName','Ankle');
%             hold on
%             plot(R.COPR(ictt1,1)*1e3,R.COPR(ictt1,2)*1e3,'o','DisplayName','COP');
%             plot(grfx',grfy','-','DisplayName','GRF');
%             grid on
%             axis equal
%             xlabel('fore-after (mm)')
%             ylabel('vertical (mm)')
%             title('Centre of pressure during stance')
            
        end
        
        
    end
    
    
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
    set(h5,'PaperPositionMode','auto')
    print(h5,[figNamePrefix '_qs_lit'],'-dpng','-r0')
    set(h6,'PaperPositionMode','auto')
    print(h6,[figNamePrefix '_P_lit'],'-dpng','-r0')
    set(h7,'PaperPositionMode','auto')
    print(h7,[figNamePrefix '_F_PF_lit'],'-dpng','-r0')
    set(h8,'PaperPositionMode','auto')
    print(h8,[figNamePrefix '_COP'],'-dpng','-r0')
end

end


