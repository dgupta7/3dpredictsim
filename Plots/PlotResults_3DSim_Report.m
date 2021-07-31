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
fsq = [scs(3)/2, scs(4)*0.6];
fhigh = [scs(3)/2, scs(4)-120];
flong = [scs(3)/2, scs(4)/4];
fhigh1 = [scs(3)*0.4, scs(4)-120];

fpos = [1,scs(4)/2+20;
        1,40;
        scs(3)/2,40;
        0,40;];

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
    Cs = CsV(inr,:);
    
    if inr==1
        hleg = figure;
        hold on
        plot(0,R.COT,'o','Color',Cs,'MarkerFaceColor',Cs,'DisplayName',LegName);
        lh=legend('-DynamicLegend','location','northwestoutside');
        title(lh,'Legend')
        title('Cost of transport')
        ylabel('COT (J/kg/m)')
    else
        figure(hleg)
        plot(0,R.COT,'o','Color',Cs,'MarkerFaceColor',Cs,'DisplayName',LegName);

        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(hleg,'PaperPositionMode','auto')
            print(hleg,[figNamePrefix '_legend'],'-dpng','-r0')
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
                lh=legend('-DynamicLegend','location','northwest');
                lh.Interpreter = 'none';
                lhPos = lh.Position;
                lhPos(1) = lhPos(1)+0.15;
                lhPos(2) = lhPos(2)+0.1;
                set(lh,'position',lhPos);
            end
        end
        
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h1,'PaperPositionMode','auto')
            print(h1,[figNamePrefix '_qs'],'-dpng','-r0')
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
                lhPos(1) = lhPos(1)+0.15;
                lhPos(2) = lhPos(2)+0.1;
                set(lh,'position',lhPos);
            end
        end
        
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h2,'PaperPositionMode','auto')
            print(h2,[figNamePrefix '_Ts'],'-dpng','-r0')
        end
    end
    
    
    %% soleus
    iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
    iGas = find(strcmp(R.colheaders.muscles,'lat_gas_r'));
    iGas2 = find(strcmp(R.colheaders.muscles,'med_gas_r'));
    if isempty(iGas)
        iGas = find(strcmp(R.colheaders.muscles,'gaslat_r'));
    end
    if isempty(iGas2)
        iGas2 = find(strcmp(R.colheaders.muscles,'gasmed_r'));
    end
        
    if makeplot.soleus
        
        if inr==1
            h3 = figure('Position',[fpos(3,:),fhigh]);
        else
            figure(h3);
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
                subplot(7,3,imu); hold on;
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
        NumTicks = 6;
        
        subplot(7,3,1); hold on;
        plot(R.a(:,iSol),'-','Color',CsV(inr,:),'DisplayName',LegName);
        title('Soleus')
        ylabel('Activity (-)')
        grid on
        L = get(gca,'XLim');
        set(gca,'XTick',linspace(L(1),L(2),NumTicks))
        
        subplot(7,3,2); hold on;
        plot(R.a(:,iGas),'-','Color',CsV(inr,:),'DisplayName',LegName);
        title('Gastrocnemius-medialis')
        grid on
        L = get(gca,'XLim');
        set(gca,'XTick',linspace(L(1),L(2),NumTicks))
        
        subplot(7,3,3); hold on;
        plot(R.a(:,iGas2),'-','Color',CsV(inr,:),'DisplayName',LegName);
        title('Gastrocnemius-lateralis')
        grid on
        L = get(gca,'XLim');
        set(gca,'XTick',linspace(L(1),L(2),NumTicks))
        
        if inr==1
            lh=legend('-DynamicLegend','location','northeast');
            lh.Interpreter = 'none';
        end
        
        
        for imu=1:3
            
            subplot(7,3,3+imu)
            plot(R.FT(:,imus(imu)),'-','Color',CsV(inr,:),'DisplayName',LegName);
            hold on
            grid on
            
            L = get(gca,'XLim');
            set(gca,'XTick',linspace(L(1),L(2),NumTicks))
            if imu==1
                ylabel('F_{normal} (N)');
            end

            subplot(7,3,6+imu)
            plot(R.MetabB.Etot(:,imus(imu)),'-','Color',CsV(inr,:)); hold on;
            hold on
            grid on
            
            L = get(gca,'XLim');
            set(gca,'XTick',linspace(L(1),L(2),NumTicks))
            if imu==1
                ylabel('P_{metabolic} (W)');
            end

            subplot(7,3,9+imu)
            plot(R.MetabB.Wdot(:,imus(imu)),'-','Color',CsV(inr,:)); hold on;
            hold on
            grid on
            L = get(gca,'XLim');
            set(gca,'XTick',linspace(L(1),L(2),NumTicks))
            if imu==1
                ylabel('P_{mech,tot} (W)');
            end
            
            if isfield(R,'vT')
               
                subplot(7,3,12+imu)
                plot(-R.FT(:,imus(imu)).*R.vT(:,imus(imu)),'-','Color',CsV(inr,:)); hold on;
                hold on
                grid on
                
                L = get(gca,'XLim');
                set(gca,'XTick',linspace(L(1),L(2),NumTicks))
                if imu==1
                    ylabel('P_{tendon} (W)');
                end
            end
       
            subplot(7,3,15+imu)
            plot(R.lMtilde(:,imus(imu)),'-','Color',CsV(inr,:)); hold on;
            hold on
            grid on
            L = get(gca,'XLim');
            set(gca,'XTick',linspace(L(1),L(2),NumTicks))
            if imu==1
                ylabel('Fibre length (-)');
            end

            subplot(7,3,18+imu)
            plot(R.Muscle.vM(:,imus(imu)),'Color',CsV(inr,:),'DisplayName',LegName);
            hold on
            grid on
            L = get(gca,'XLim');
            set(gca,'XTick',linspace(L(1),L(2),NumTicks))
            if imu==1
                ylabel('Fibre velocity (s^{-1})')
            end

            xlabel('Gait cycle (%)','Fontsize',label_fontsize);
        end
        lhPos = lh.Position;
        lhPos(1) = lhPos(1)+0.1;
%         lhPos(2) = lhPos(2)+0.08;
        set(lh,'position',lhPos);
        
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h3,'PaperPositionMode','auto')
            print(h3,[figNamePrefix '_calf'],'-dpng','-r0')
        end
    end
    
    
    
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
            L = get(gca,'XLim');
            set(gca,'XTick',linspace(L(1),L(2),NumTicks))
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
        
%         figure
%         plot(R.GRFs(:,2)+R.GRFs(:,5),'-')
%         hold on
%         plot(Dat.(type).gc.GRF.Fmean(:,2)+Dat.(type).gc.GRF.Fmean([51:end,1:50],2),'-k')
    
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h4,'PaperPositionMode','auto')
            print(h4,[figNamePrefix '_GRF'],'-dpng','-r0')
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
            
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h5,'PaperPositionMode','auto')
            print(h5,[figNamePrefix '_qs_lit'],'-dpng','-r0')
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
        
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h8,'PaperPositionMode','auto')
            print(h8,[figNamePrefix '_COP'],'-dpng','-r0')
        end
        
    end
    
    
    %% all
    if makeplot.allQsTs
        if inr==1
            h9 = figure('Position',[fpos(3,:),fhigh1]);
            h10 = figure('Position',[fpos(4,:),fhigh1]);
        end
        
        if has_no_tmt && has_no_mtj
            idx_Qs = [1,2,3,10,11,12,14,16,18,20,21,22,23,27,28,29,31];
        else
            idx_Qs = [1,2,3,10,11,12,14,16,18,20,22,23,24,25,29,30,31,33];
        end
        idx_title = [1,2,3,10,11,12,14,16,18,20,24,25,26,27,31,32,33,35];
        joints_ref = {'pelvis_tilt','pelvis_list','pelvis_rotation',...
                'hip_flexion','hip_adduction','hip_rotation',...
                'knee_angle','ankle_angle','subtalar_angle','mtj_angle','mtp_angle',...
                'lumbar_extension','lumbar_bending','lumbar_rotation',...
                'arm_flex','arm_add','arm_rot','elbow_flex'};
        
        joints_tit = {'Pelvis tilt','Pelvis list','Pelvis rotation','Pelvis tx',...
                'Pelvis ty','Pelvis tz','Hip flexion L','Hip adduction L',...
                'Hip rotation L','Hip flexion R','Hip adduction R','Hip rotation R',...
                'Knee L','Knee R','Ankle L','Ankle R','Subtalar L','Subtalar R',...
                'Midtarsal L','Midtarsal R','Tmt L','Tmt R',...
                'Mtp L','Mtp R',...
                'Lumbar extension','Lumbar bending','Lumbar rotation',...
                'Arm flexion L','Arm adduction L','Arm rotation L',...
                'Arm flexion R','Arm adduction R','Arm rotation R',...
                'Elbow flexion L','Elbow flexion R'};
    
        figure(h9)
        
        j = 0;
        label_fontsize  = 12;
        line_linewidth  = 0.5;
        for i = 1:length(idx_title)
            subplot(6,3,i)
            x = 1:(100-1)/(size(R.Qs,1)-1):100;
            % Experimental data
            if  inr == 1 && md
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
                    fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k','DisplayName',['MoCap ' refName]);
                    alpha(.25);
                end
            end

            % Simulation results
            x = 1:(100-1)/(size(R.Qs,1)-1):100;
            hold on;
            axis tight
            xlim([0,100]);
            if (has_no_tmt && strcmp(joints_tit{idx_title(i)},'Tarsometatarsal R')) || ...
                    (has_no_mtj && strcmp(joints_tit{idx_title(i)},'Midtarsal R'))
                % skip this plot
            else
                j=j+1;
                plot(x,R.Qs(:,idx_Qs(j)),'color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
            end

            % Plot settings
            if inr==1
                set(gca,'Fontsize',label_fontsize);
                title(joints_tit{idx_title(i)},'Fontsize',label_fontsize);
                % Y-axis
                if i == 1 || i == 4 || i == 7 || i == 10 || i == 13 || i == 16
                    ylabel('Angle (°)','Fontsize',label_fontsize);
                end
                % X-axis
                L = get(gca,'XLim');
                NumTicks = 3;
                set(gca,'XTick',linspace(L(1),L(2),NumTicks))
                if i > 15
                    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
                end
            end
            if inr==1 && i==3
                lh=legend('-DynamicLegend','location','west');
                lh.Interpreter = 'none';
                lhPos = lh.Position;
                lhPos(2) = lhPos(2)+0.085;
                set(lh,'position',lhPos);
            end
        end

        
        figure(h10)
        j = 0;
        label_fontsize  = 12;
        line_linewidth  = 0.5;
        for i = 1:length(idx_title)
            subplot(6,3,i)
            x = 1:(100-1)/(size(R.Qs,1)-1):100;
            % Experimental data
            if  inr == 1 && md
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
                    fill([x fliplr(x)],[meanPlusSTD fliplr(meanMinusSTD)],'k','DisplayName',['MoCap ' refName]);
                    alpha(.25);
                end
            end

            % Simulation results
            x = 1:(100-1)/(size(R.Qs,1)-1):100;
            hold on;
            axis tight
            xlim([0,100]);
            if (has_no_tmt && strcmp(joints_tit{idx_title(i)},'Tarsometatarsal R')) || ...
                    (has_no_mtj && strcmp(joints_tit{idx_title(i)},'Midtarsal R'))
                % skip this plot
            else
                j=j+1;
                plot(x,R.Tid(:,idx_Qs(j)),'color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
            end
            % Plot settings
            if inr==1
                % Plot settings
                set(gca,'Fontsize',label_fontsize);
                title(joints_tit{idx_title(i)},'Fontsize',label_fontsize);
                % Y-axis
                if i == 1 || i == 4 || i == 7 || i == 10 || i == 13 || i == 16
                    ylabel('Torque (Nm)','Fontsize',label_fontsize);
                end
                % X-axis
                L = get(gca,'XLim');
                NumTicks = 3;
                set(gca,'XTick',linspace(L(1),L(2),NumTicks))
                if i > 15
                    xlabel('Gait cycle (%)','Fontsize',label_fontsize);
                end
            end
            if inr==1 && i==3
                lh=legend('-DynamicLegend','location','west');
                lh.Interpreter = 'none';
                lhPos = lh.Position;
                lhPos(2) = lhPos(2)+0.085;
                set(lh,'position',lhPos);
            end
        end 
        
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h9,'PaperPositionMode','auto')
            print(h9,[figNamePrefix '_qs_all'],'-dpng','-r0')
            set(h10,'PaperPositionMode','auto')
            print(h10,[figNamePrefix '_Ts_all'],'-dpng','-r0')
        end
        
    end
        
    
    %%
    x = 1:(100-1)/(size(R.Qs,1)-1):100;
    if makeplot.k_mtj_lin
        if inr==1
            h11 = figure('Position',[fpos(4,:),fwide]);
        end
        
        imtj = find(strcmp(R.colheaders.joints,'mtj_angle_r'));
        iankle = strcmp(R.colheaders.joints,'ankle_angle_r');
        isubt = strcmp(R.colheaders.joints,'subtalar_angle_r');
        imtp = find(strcmp(R.colheaders.joints,'mtp_angle_r'));
        
    
        if ~isempty(imtj) && (R.S.MT_li_nonl == 0 || strcmp(R.S.mtj_stiffness,'signed_lin'))
        
            W_ankle = zeros(length(x),1);
            W_subt = zeros(length(x),1);
            W_mtp = zeros(length(x),1);
            W_mtj = zeros(length(x),1);
            
            for ii=2:length(x)
               W_ankle(ii) = trapz(R.Qs(1:ii,iankle)*pi/180,R.Tid(1:ii,iankle));
               W_subt(ii) = trapz(R.Qs(1:ii,isubt)*pi/180,R.Tid(1:ii,isubt));
               W_mtp(ii) = trapz(R.Qs(1:ii,imtp)*pi/180,R.Tid(1:ii,imtp));
               W_mtj(ii) = trapz(R.Qs(1:ii,imtj)*pi/180,R.Tid(1:ii,imtj));
            end
            W_mtj(2:end) = (W_mtj(2:end)-W_mtj(2))/R.body_mass;
            W_ankle(2:end) = (W_ankle(2:end)-W_ankle(2))/R.body_mass;
            W_subt(2:end) = (W_subt(2:end)-W_subt(2))/R.body_mass;
            W_mtp(2:end) = (W_mtp(2:end)-W_mtp(2))/R.body_mass;
            W_tot = W_mtj + W_ankle + W_subt + W_mtp;
            
            figure(h11)
            subplot(2,4,1)
            hold on
            grid on
            plot(R.S.kMT_li,W_tot(end),'o','Color',Cs,'MarkerFaceColor',Cs)
            ylabel('Work (J/kg)')
            title('W_{mech} foot (1 GC)')

            dist_trav = R.Qs(end,strcmp(R.colheaders.joints,'pelvis_tx')) - ...
                R.Qs(1,strcmp(R.colheaders.joints,'pelvis_tx'));

            subplot(2,4,5)
            hold on
            grid on
            plot(R.S.kMT_li,R.COT*dist_trav,'o','Color',Cs,'MarkerFaceColor',Cs)
            xlabel('mtj stiffness (Nm/rad)')
            ylabel('Metab. Energy (J/kg)')
            title('E_{metab} foot (1 GC)')

            subplot(2,4,2)
            hold on
            grid on
            plot(R.S.kMT_li,W_tot(end)/dist_trav,'o','Color',Cs,'MarkerFaceColor',Cs)
            ylabel('Work (J/(kg m))')
            title({'W_{mech} foot (1 m)'})

            subplot(2,4,6)
            hold on
            grid on
            plot(R.S.kMT_li,R.COT,'o','Color',Cs,'MarkerFaceColor',Cs)
            xlabel('mtj stiffness (Nm/rad)')
            ylabel('COT (J/(kg m)')
            title('Cost of transport')


            subplot(2,4,3)
            hold on
            grid on
            plot(R.S.kMT_li,max(R.Qs(:,imtj)),'o','Color',Cs,'MarkerFaceColor',Cs)
            ylabel('angle (°)')
            title({'Max mtj angle'})
        
            subplot(2,4,7)
            hold on
            grid on
            plot(R.S.kMT_li,min(R.Qs(:,imtj)),'o','Color',Cs,'MarkerFaceColor',Cs)
            xlabel('mtj stiffness (Nm/rad)')
            ylabel('angle (°)')
            title({'Min mtj angle'})
            
            subplot(2,4,4)
            hold on
            grid on
            plot(R.S.kMT_li,max(R.Qs(:,imtp)),'o','Color',Cs,'MarkerFaceColor',Cs)
            ylabel('angle (°)')
            title({'Max mtp angle'})
        
            subplot(2,4,8)
            hold on
            grid on
            plot(R.S.kMT_li,min(R.Qs(:,imtp)),'o','Color',Cs,'MarkerFaceColor',Cs)
            xlabel('mtj stiffness (Nm/rad)')
            ylabel('angle (°)')
            title({'Min mtp angle'})
            
        end
        
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h11,'PaperPositionMode','auto')
            print(h11,[figNamePrefix '_k_mtj_lin'],'-dpng','-r0')
        end
        
    end
       
    %% Windlass
    if makeplot.windlass
        if inr==1
            h12 = figure('Position',[fpos(4,:),fwide]);
        end
        
        imtj = find(strcmp(R.colheaders.joints,'mtj_angle_r'));
        iankle = strcmp(R.colheaders.joints,'ankle_angle_r');
        isubt = strcmp(R.colheaders.joints,'subtalar_angle_r');
        imtp = find(strcmp(R.colheaders.joints,'mtp_angle_r'));
        
        
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h12,'PaperPositionMode','auto')
            print(h12,[figNamePrefix '_WL'],'-dpng','-r0')
        end
    end
    
    %% foot power
    if makeplot.power || makeplot.work
        
        x = 1:(100-1)/(size(R.Qs,1)-1):100;
        istance = 1:1:ceil(R.Event.Stance)+10;
        xst = linspace(1,110,length(istance));
        
        imtj = find(strcmp(R.colheaders.joints,'mtj_angle_r'));
        iknee = strcmp(R.colheaders.joints,'knee_angle_r');
        iankle = strcmp(R.colheaders.joints,'ankle_angle_r');
        isubt = strcmp(R.colheaders.joints,'subtalar_angle_r');
        imtp = find(strcmp(R.colheaders.joints,'mtp_angle_r'));
        
        P_knee = R.Qdots(:,iknee)*pi/180.*R.Tid(:,iknee)/R.body_mass;
        P_ankle = R.Qdots(:,iankle)*pi/180.*R.Tid(:,iankle)/R.body_mass;
        P_subt = R.Qdots(:,isubt)*pi/180.*R.Tid(:,isubt)/R.body_mass;
        P_mtp = R.Qdots(:,imtp)*pi/180.*R.Tid(:,imtp)/R.body_mass;
        
        W_knee = zeros(size(P_knee));
        W_ankle = zeros(size(P_ankle));
        W_subt = zeros(size(P_subt));
        W_mtp = zeros(size(P_mtp));
        
        for iw=2:length(xst)
           W_knee(iw) = trapz(R.t(1:iw),P_knee(1:iw));
           W_ankle(iw) = trapz(R.t(1:iw),P_ankle(1:iw));
           W_subt(iw) = trapz(R.t(1:iw),P_subt(1:iw));
           W_mtp(iw) = trapz(R.t(1:iw),P_mtp(1:iw));
        end
        
        if ~isempty(imtj)
            P_mtj = R.Qdots(:,imtj)*pi/180.*R.Tid(:,imtj)/R.body_mass;
            W_mtj = zeros(size(P_mtj));
            
            l_PF = R.windlass.l_PF;
            F_PF = R.windlass.F_PF;
            qdot_mtj = R.Qdots(:,imtj)*pi/180;
            M_li = R.windlass.M_li;
            M_mtp = R.windlass.M_mtp;
            
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
            P_mtp_WL = R.Qdots(:,imtp)*pi/180.*M_mtp/R.body_mass;
        
            W_PF = zeros(size(P_PF));
            W_li = zeros(size(P_li));
            W_mtp_WL = zeros(size(P_mtp_WL));
            
            for iw=2:length(xst)
                W_PF(iw) = trapz(R.t(1:iw),P_PF(1:iw));
                W_li(iw) = trapz(R.t(1:iw),P_li(1:iw));
                W_mtp_WL(iw) = trapz(R.t(1:iw),P_mtp_WL(1:iw));
                W_mtj(iw) = trapz(R.t(1:iw),P_mtj(1:iw));
            end
            
        else
            P_mtj = zeros(size(P_mtp));
            W_mtj = zeros(size(W_mtp));
        end
        
        P_joints = P_mtj + P_ankle + P_subt + P_mtp;
        W_joints = W_mtj + W_ankle + W_subt + W_mtp;

        if isfield(R,'vT')
           P_T_Sol = -R.FT(:,iSol).*R.vT(:,iSol)/R.body_mass;
           P_T_Gas = -R.FT(:,iGas).*R.vT(:,iGas)/R.body_mass;
           P_T_Gas2 = -R.FT(:,iGas2).*R.vT(:,iGas2)/R.body_mass;
           
           W_T_Sol = zeros(size(P_T_Sol));
           W_T_Gas = zeros(size(P_T_Sol));
           W_T_Gas2 = zeros(size(P_T_Sol));
           
           for iw=2:length(xst)
                W_T_Sol(iw) = trapz(R.t(1:iw),P_T_Sol(1:iw));
                W_T_Gas(iw) = trapz(R.t(1:iw),P_T_Gas(1:iw));
                W_T_Gas2(iw) = trapz(R.t(1:iw),P_T_Gas2(1:iw));
           end
           
        end
    
        P_M_Sol = R.MetabB.Wdot(:,iSol)/R.body_mass;
        P_M_Gas = R.MetabB.Wdot(:,iGas)/R.body_mass;
        P_M_Gas2 = R.MetabB.Wdot(:,iGas2)/R.body_mass;
        
        W_M_Sol = zeros(size(P_M_Sol));
        W_M_Gas = zeros(size(P_M_Sol));
        W_M_Gas2 = zeros(size(P_M_Sol));
       
        for iw=2:length(xst)
            W_M_Sol(iw) = trapz(R.t(1:iw),P_M_Sol(1:iw));
            W_M_Gas(iw) = trapz(R.t(1:iw),P_M_Gas(1:iw));
            W_M_Gas2(iw) = trapz(R.t(1:iw),P_M_Gas2(1:iw));
        end
           
        P_HC_heel = R.P_mech_contact.vertical.calcn.r/R.body_mass;
        P_HC_ball = R.P_mech_contact.vertical.metatarsi.r/R.body_mass;
        P_HC_toes = R.P_mech_contact.vertical.toes.r/R.body_mass;
        P_HC = P_HC_heel + P_HC_ball + P_HC_toes;
        
        W_HC_heel = zeros(size(P_HC));
        W_HC_ball = zeros(size(P_HC));
        W_HC_toes = zeros(size(P_HC));
        
        for iw=2:length(xst)
            W_HC_heel(iw) = trapz(R.t(1:iw),R.P_mech_contact.vertical.calcn.r(1:iw))/R.body_mass;
            W_HC_ball(iw) = trapz(R.t(1:iw),R.P_mech_contact.vertical.metatarsi.r(1:iw))/R.body_mass;
            W_HC_toes(iw) = trapz(R.t(1:iw),R.P_mech_contact.vertical.toes.r(1:iw))/R.body_mass;
        end
        
        W_HC = W_HC_heel + W_HC_ball + W_HC_toes;
        P_tot = P_HC + P_joints;
        W_tot = W_HC + W_joints;
    end
    


%%
    if makeplot.power
        if inr==1
            h14 = figure('Position',[fpos(4,:),fsq]);
        end
        
        figure(h14)
        
        nh = 3;
        nw = 4;
        
        pos_P_all = {P_tot, P_joints, P_HC, 'none',...
                   P_HC_heel, P_HC_ball, P_HC_toes, 'none',...
                   P_ankle, P_subt, P_mtj, P_mtp};
               
        titles_P_all = {{'Ankle-foot','\fontsize{10}\rm joints + pads'},...
            {'Joints','\fontsize{10}\rm ankle + subt + mtj + mtp'},...
            {'Pads','\fontsize{10}\rm hindfoot + forefoot + toes'},...
            'none', 'Hindfoot','Forefoot','Toes','none',...
            'Ankle','Subt','Mtj','Mtp'};
        
        for i_P=1:numel(pos_P_all)
            if ~strcmp(pos_P_all{i_P},'none')
                subplot(nh,nw,i_P)
                hold on
                grid on
                plot(xst,pos_P_all{i_P}(istance),'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
                title(titles_P_all{i_P},'Fontsize',label_fontsize);
                axis tight
                yl = get(gca, 'ylim');
                ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)])
                xlim([0,110,])
                
                if mod(i_P-1,nw)==0
                    ylabel('P_{mech} (W/kg)','Fontsize',label_fontsize);
                end
                if i_P > nh*nw-nw
                    xlabel('Stance phase (%)','Fontsize',label_fontsize);
                end
            end
        end
        
       if isfield(R,'GRFs_separate') && ~isempty(R.GRFs_separate)
            subplot(nh,nw,8)
            hold on
            grid on
            p1=plot(xst,R.GRFs_separate(istance,2),'-.','Color',Cs,'DisplayName','calcaneus');
            p2=plot(xst,R.GRFs_separate(istance,2+3),'--','Color',Cs,'DisplayName','forefoot');
            p3=plot(xst,R.GRFs_separate(istance,2+6),':','Color',Cs,'DisplayName','toes');
            xlabel('% stride');
            ylabel('% body weight')
            title('Vertical GRF')
            lh=legend([p1,p2,p3],'location','best');
            lh.Interpreter = 'none';
        end
        
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h14,'PaperPositionMode','auto')
            print(h14,[figNamePrefix '_P1'],'-dpng','-r0')
        end
        
        %%
        if inr==1
            h15 = figure('Position',[fpos(3,:),fsq]);
        end
        
        figure(h15)
        
        nh = 3;
        nw = 4;
        
        pos_P_all = {P_knee, P_ankle, P_subt, P_knee+P_ankle+P_subt,...
            P_T_Sol, P_T_Gas, P_T_Gas2,P_T_Sol+P_T_Gas+P_T_Gas2,...
            P_M_Sol, P_M_Gas, P_M_Gas2,P_M_Sol+P_M_Gas+P_M_Gas2};
        
        
               
        titles_P_all = {'Knee','Ankle','Subt',{'Sum','\fontsize{10}\rm knee + ankle + subt'},...
            'Soleus tendon','Gas-med tendon','Gas-lat tendon',{'Sum Tendon','\fontsize{10}\rm sol + gas-med + gas-lat'},...
            'Soleus fibres','Gas-med fibres','Gas-lat fibres',{'Sum fibres','\fontsize{10}\rm sol + gas-med + gas-lat'}};
        
        for i_P=1:numel(pos_P_all)
            if ~strcmp(pos_P_all{i_P},'none')
                subplot(nh,nw,i_P)
                hold on
                grid on
                plot(xst,pos_P_all{i_P}(istance),'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
                title(titles_P_all{i_P},'Fontsize',label_fontsize);
                axis tight
                yl = get(gca, 'ylim');
                ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)])
                xlim([0,110,])
                
                if mod(i_P-1,nw)==0
                    ylabel('P_{mech} (W/kg)','Fontsize',label_fontsize);
                end
                if i_P > nh*nw-nw
                    xlabel('Stance phase (%)','Fontsize',label_fontsize);
                end
            end
        end
        
        if isfield(R,'GRFs_separate') && ~isempty(R.GRFs_separate)
            subplot(nh,nw,4)
            hold on
            grid on
            plot(xst,R.GRFs_separate(istance,2),'Color',Cs);
            ylabel('% body weight','Fontsize',label_fontsize);
            title('Vertical GRF Heel')
            xlim([0,110,])

            subplot(nh,nw,8)
            hold on
            grid on
            plot(xst,R.GRFs_separate(istance,2+3),'Color',Cs);
            ylabel('% body weight','Fontsize',label_fontsize);
            title('Vertical GRF Forefoot')
            xlim([0,110,])
            
            subplot(nh,nw,12)
            hold on
            grid on
            plot(xst,R.GRFs_separate(istance,2+6),'Color',Cs);
            xlabel('Stance phase (%)','Fontsize',label_fontsize);
            ylabel('% body weight','Fontsize',label_fontsize);
            title('Vertical GRF Toes')
            xlim([0,110,])
        end
        
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h15,'PaperPositionMode','auto')
            print(h15,[figNamePrefix '_P_ankle'],'-dpng','-r0')
        end
        
%%
        if inr==1
            h16 = figure('Position',[fpos(3,:),fsq]);
        end
        
        figure(h16)
        
        nh = 3;
        nw = 4;
        
        if ~isempty(imtj)
            pos_P_all = {P_mtj, P_mtp, P_mtj+P_mtp, 'none',...
                P_mtj-P_li, P_mtp_WL, P_PF, 'none',...
                P_li, P_mtp-P_mtp_WL, 'none', 'none'};
        else
            pos_P_all = {P_mtj, P_mtp, P_mtj+P_mtp, 'none',...
                'none', 'none', 'none', 'none',...
                'none', 'none', 'none', 'none'};
        end
        
               
        titles_P_all = {'Mtj','Mtp',{'Sum','\fontsize{10}\rm mtj + mtj'},'none'...
            'Mtj PF','Mtp PF','PF','none', 'Mtj non-PF','Mtp non-PF','none','none'};
        
        for i_P=1:numel(pos_P_all)
            if ~strcmp(pos_P_all{i_P},'none')
                subplot(nh,nw,i_P)
                hold on
                grid on
                plot(xst,pos_P_all{i_P}(istance),'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
                title(titles_P_all{i_P},'Fontsize',label_fontsize);
                axis tight
                yl = get(gca, 'ylim');
                ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)])
                xlim([0,110,])
                
                if mod(i_P-1,nw)==0
                    ylabel('P_{mech} (W/kg)','Fontsize',label_fontsize);
                end
                if i_P > nh*nw-nw
                    xlabel('Stance phase (%)','Fontsize',label_fontsize);
                end
            end
        end
        
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h16,'PaperPositionMode','auto')
            print(h16,[figNamePrefix '_P_foot'],'-dpng','-r0')
        end
        
        if makeplot.power_T
            P_dist_hindfoot = P_mtp + P_mtj + P_HC;
            P_dist_forefoot = P_mtp + (R.P_mech_contact.vertical.metatarsi.r+R.P_mech_contact.vertical.toes.r)/R.body_mass;
            P_dist_hallux = R.P_mech_contact.vertical.toes.r/R.body_mass;

            figure
            subplot(8,1,1:5)
            hold on
            grid on
            plot(xst,P_dist_hallux(istance),'-','Color',[0, 0.5, 0],'linewidth',line_linewidth,'DisplayName','Hallux');
            plot(xst,P_dist_forefoot(istance),'-','Color',[0.7, 0.1, 0.1],'linewidth',line_linewidth,'DisplayName','Forefoot');
            plot(xst,P_dist_hindfoot(istance),'-','Color',[0, 0.4470, 0.7410],'linewidth',line_linewidth,'DisplayName','Hindfoot');
            plot(xst,P_tot(istance),'-','Color','k','linewidth',line_linewidth,'DisplayName','Shank');
            ylabel('P_{mech} (W/kg)')
            xlabel('Gait cycle (%)','Fontsize',label_fontsize);
    %         title('Total foot')
            xlim([0,110,])
            lg = legend('Location','southeast');
            title(lg,'Distal to ...')

            W_tot = trapz(R.t,P_tot);
            W_dist_hindfoot = trapz(R.t,P_dist_hindfoot);
            W_dist_forefoot = trapz(R.t,P_dist_forefoot);
            W_dist_hallux = trapz(R.t,P_dist_hallux);

            W_s = [W_tot,W_dist_hindfoot,W_dist_forefoot,W_dist_hallux];
            c = categorical({'Shank','Hindfoot','Forefoot','Hallux'});
            c = reordercats(c,[2;1;3;4]);

            subplot(8,1,7:8)
            br=bar(c,W_s);
            grid on
            br.FaceColor = 'flat';
            br.CData(1,:) = [0, 0.5, 0];
            br.CData(2,:) = [0.7, 0.1, 0.1];
            br.CData(3,:) = [0, 0.4470, 0.7410];
            br.CData(4,:) = [0, 0, 0];
            ylabel('W_{mech,net} (J/kg)')

        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%
    %%
    if makeplot.work
        if inr==1
            h17 = figure('Position',[fpos(4,:),fsq]);
        end
        
        figure(h17)
        
        nh = 3;
        nw = 4;
        
        pos_P_all = {W_tot, W_joints, W_HC, 'none',...
                   W_HC_heel, W_HC_ball, W_HC_toes, 'none',...
                   W_ankle, W_subt, W_mtj, W_mtp};
               
        titles_P_all = {{'Ankle-foot','\fontsize{10}\rm joints + pads'},...
            {'Joints','\fontsize{10}\rm ankle + subt + mtj + mtp'},...
            {'Pads','\fontsize{10}\rm hindfoot + forefoot + toes'},...
            'none', 'Hindfoot','Forefoot','Toes','none',...
            'Ankle','Subt','Mtj','Mtp'};
        
        for i_P=1:numel(pos_P_all)
            if ~strcmp(pos_P_all{i_P},'none')
                subplot(nh,nw,i_P)
                hold on
                grid on
                plot(xst,pos_P_all{i_P}(istance),'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
                title(titles_P_all{i_P},'Fontsize',label_fontsize);
                axis tight
                yl = get(gca, 'ylim');
                ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)])
                xlim([0,110,])
                
                if mod(i_P-1,nw)==0
                    ylabel('W_{mech} (J/kg)','Fontsize',label_fontsize);
                end
                if i_P > nh*nw-nw
                    xlabel('Stance phase (%)','Fontsize',label_fontsize);
                end
            end
        end
        
       if isfield(R,'GRFs_separate') && ~isempty(R.GRFs_separate)
            subplot(nh,nw,8)
            hold on
            grid on
            p1=plot(xst,R.GRFs_separate(istance,2),'-.','Color',Cs,'DisplayName','calcaneus');
            p2=plot(xst,R.GRFs_separate(istance,2+3),'--','Color',Cs,'DisplayName','forefoot');
            p3=plot(xst,R.GRFs_separate(istance,2+6),':','Color',Cs,'DisplayName','toes');
            ylabel('% body weight')
            title('Vertical GRF')
            lh=legend([p1,p2,p3],'location','best');
            lh.Interpreter = 'none';
        end
        
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h17,'PaperPositionMode','auto')
            print(h17,[figNamePrefix '_W1'],'-dpng','-r0')
        end
        
        %%
        if inr==1
            h18 = figure('Position',[fpos(3,:),fsq]);
        end
        
        figure(h18)
        
        nh = 3;
        nw = 4;
        
        pos_P_all = {W_knee, W_ankle, W_subt, W_knee+W_ankle+W_subt,...
            W_T_Sol, W_T_Gas, W_T_Gas2,W_T_Sol+W_T_Gas+W_T_Gas2,...
            W_M_Sol, W_M_Gas, W_M_Gas2,W_M_Sol+W_M_Gas+W_M_Gas2};
        
        
               
        titles_P_all = {'Knee','Ankle','Subt',{'Sum','\fontsize{10}\rm knee + ankle + subt'},...
            'Soleus tendon','Gas-med tendon','Gas-lat tendon',{'Sum Tendon','\fontsize{10}\rm sol + gas-med + gas-lat'},...
            'Soleus fibres','Gas-med fibres','Gas-lat fibres',{'Sum fibres','\fontsize{10}\rm sol + gas-med + gas-lat'}};
        
        for i_P=1:numel(pos_P_all)
            if ~strcmp(pos_P_all{i_P},'none')
                subplot(nh,nw,i_P)
                hold on
                grid on
                plot(xst,pos_P_all{i_P}(istance),'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
                title(titles_P_all{i_P},'Fontsize',label_fontsize);
                axis tight
                yl = get(gca, 'ylim');
                ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)])
                xlim([0,110,])
                
                if mod(i_P-1,nw)==0
                    ylabel('W_{mech} (J/kg)','Fontsize',label_fontsize);
                end
                if i_P > nh*nw-nw
                    xlabel('Stance phase (%)','Fontsize',label_fontsize);
                end
            end
        end
        
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h18,'PaperPositionMode','auto')
            print(h18,[figNamePrefix '_W_ankle'],'-dpng','-r0')
        end
        
%%
        if inr==1
            h19 = figure('Position',[fpos(3,:),fsq]);
        end
        
        figure(h19)
        
        nh = 3;
        nw = 4;
        
        if ~isempty(imtj)
            pos_P_all = {W_mtj, W_mtp, W_mtj+W_mtp, 'none',...
                W_mtj-W_li, W_mtp_WL, W_PF, 'none',...
                W_li, W_mtp-W_mtp_WL, 'none', 'none'};
        else
            pos_P_all = {W_mtj, W_mtp, W_mtj+W_mtp, 'none',...
                'none', 'none', 'none', 'none',...
                'none', 'none', 'none', 'none'};
        end
        
               
        titles_P_all = {'Mtj','Mtp',{'Sum','\fontsize{10}\rm mtj + mtj'},'none'...
            'Mtj PF','Mtp PF','PF','none', 'Mtj non-PF','Mtp non-PF','none','none'};
        
        for i_P=1:numel(pos_P_all)
            if ~strcmp(pos_P_all{i_P},'none')
                subplot(nh,nw,i_P)
                hold on
                grid on
                plot(xst,pos_P_all{i_P}(istance),'Color',Cs,'linewidth',line_linewidth,'DisplayName',LegName);
                title(titles_P_all{i_P},'Fontsize',label_fontsize);
                axis tight
                yl = get(gca, 'ylim');
                ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)])
                xlim([0,110,])
                
                if mod(i_P-1,nw)==0
                    ylabel('W_{mech} (J/kg)','Fontsize',label_fontsize);
                end
                if i_P > nh*nw-nw
                    xlabel('Stance phase (%)','Fontsize',label_fontsize);
                end
            end
        end
        
        if isfield(R,'GRFs_separate') && ~isempty(R.GRFs_separate)
            subplot(nh,nw,4)
            hold on
            grid on
            plot(xst,R.GRFs_separate(istance,2),'Color',Cs);
            ylabel('% body weight','Fontsize',label_fontsize);
            title('Vertical GRF Heel')
            axis tight
            yl = get(gca, 'ylim');
            ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)])
            xlim([0,110,])

            subplot(nh,nw,8)
            hold on
            grid on
            plot(xst,R.GRFs_separate(istance,2+3),'Color',Cs);
            ylabel('% body weight','Fontsize',label_fontsize);
            title('Vertical GRF Forefoot')
            axis tight
            yl = get(gca, 'ylim');
            ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)])
            xlim([0,110,])
            
            subplot(nh,nw,12)
            hold on
            grid on
            plot(xst,R.GRFs_separate(istance,2+6),'Color',Cs);
            xlabel('Stance phase (%)','Fontsize',label_fontsize);
            ylabel('% body weight','Fontsize',label_fontsize);
            title('Vertical GRF Toes')
            axis tight
            yl = get(gca, 'ylim');
            ylim([yl(1)-0.1*norm(yl),yl(2)+0.1*norm(yl)])
            xlim([0,110,])
        end
        
        if inr==nr && ~strcmp(figNamePrefix,'none')
            set(h19,'PaperPositionMode','auto')
            print(h19,[figNamePrefix '_W_foot'],'-dpng','-r0')
        end
        
    end
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end



end


