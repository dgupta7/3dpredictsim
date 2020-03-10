function [] = PlotResults_3DSim(ResultsFile,Cs,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if exist(ResultsFile,'file')
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
        tab6 = h.Parent.Children(1).Children(1).Children(6);
        tab7 = h.Parent.Children(1).Children(1).Children(7);
        tab8 = h.Parent.Children(1).Children(1).Children(8);
    else
        h = figure();
        set(h,'Position',get(0,'ScreenSize'));
        hTabGroup = uitabgroup;
        tab1 = uitab(hTabGroup, 'Title', 'MainInfo');
        tab2 = uitab(hTabGroup, 'Title', 'ExoT');
        tab3 = uitab(hTabGroup, 'Title', 'CalfM');
        tab4 = uitab(hTabGroup, 'Title', 'Kinematics - Sag');
        tab5 = uitab(hTabGroup, 'Title', 'Kinetics - Sag');
        tab6 = uitab(hTabGroup, 'Title', 'Kinematics - Front');
        tab7 = uitab(hTabGroup, 'Title', 'Kinetics - Front');
        tab8 = uitab(hTabGroup, 'Title', 'Ground reaction force');
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
%     subplot(2,2,3);  hold on;
%     plot(R.T_exo(:,1),'-','Color',Cs);
%     ylabel('Exo Moment [Nm]');  xlabel('% stride');
%     title('Left');
     subplot(2,2,3);  hold on;
     plot(R.Sopt.ExoScale,R.dt_exoShift,'o','Color',Cs,'MarkerFaceColor',Cs);
     ylabel('Time shift exo torque[s]');  xlabel('Level assistance');
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
    
    %% Plot joint kinematics sagital plane
    
    axes('parent', tab4);
    
    % headers
    header =  R.colheaders.joints;
    
    % Plot pelvis kinematics
    nr = 4; nc= 3;
    for i=1:6
        subplot(nr,nc,i);
        plot(R.Qs(:,i),'-','Color',Cs); hold on;
        title(header{i});
    end
    
    % plot sagital plane kinematics
    subplot(nr,nc,7);
    plot(R.Qs(:,strcmp(header,'hip_flexion_r')),'-','Color',Cs); hold on;
    title('hip flexion');
    subplot(nr,nc,8);
    plot(R.Qs(:,strcmp(header,'knee_angle_r')),'-','Color',Cs); hold on;
    title('knee flexion');
    subplot(nr,nc,9);
    plot(R.Qs(:,strcmp(header,'ankle_angle_r')),'-','Color',Cs); hold on;
    title('ankle flexion');
    subplot(nr,nc,10);
    plot(R.Qs(:,strcmp(header,'mtp_angle_r')),'-','Color',Cs); hold on;
    title('mtp angle');
    subplot(nr,nc,11);
    plot(R.Qs(:,strcmp(header,'lumbar_extension')),'-','Color',Cs); hold on;
    title('lumbar extension');
    
    for i=1:nc
        subplot(nr,nc,nr*nc-i+1)
        xlabel('% gait cycle');
    end
    for i=1:nr
        subplot(nr,nc,(i-1)*nc+1)
        ylabel('Joint kinematics [deg or m]');
    end
    
    
    %% Plot joint kinetics sagital plane
    
    axes('parent', tab5);
    
    % headers
    header =  R.colheaders.joints;
    
    % Plot pelvis kinematics
    nr = 4; nc= 3;
    for i=1:6
        subplot(nr,nc,i);
        plot(R.Tid(:,i),'-','Color',Cs); hold on;
        title(header{i});
    end
    
    % plot sagital plane kinematics
    subplot(nr,nc,7);
    plot(R.Tid(:,strcmp(header,'hip_flexion_r')),'-','Color',Cs); hold on;
    title('hip flexion');
    subplot(nr,nc,8);
    plot(R.Tid(:,strcmp(header,'knee_angle_r')),'-','Color',Cs); hold on;
    title('knee flexion');
    subplot(nr,nc,9);
    plot(R.Tid(:,strcmp(header,'ankle_angle_r')),'-','Color',Cs); hold on;
    title('ankle flexion');
    subplot(nr,nc,10);
    plot(R.Tid(:,strcmp(header,'mtp_angle_r')),'-','Color',Cs); hold on;
    title('mtp angle');
    subplot(nr,nc,11);
    plot(R.Tid(:,strcmp(header,'lumbar_extension')),'-','Color',Cs); hold on;
    title('lumbar extension');
    
    for i=1:nc
        subplot(nr,nc,nr*nc-i+1)
        xlabel('% gait cycle');
    end
    for i=1:nr
        subplot(nr,nc,(i-1)*nc+1)
        ylabel('Joint kinetics [Nm or N]');
    end
    
    
    %% Plot other kinematics
    
    axes('parent', tab6);
    
    nr = 3; nc = 3;
    subplot(nr,nc,1);
    plot(R.Qs(:,strcmp(header,'hip_flexion_r')),'-','Color',Cs); hold on;
    title('hip flexion');
    subplot(nr,nc,2);
    plot(R.Qs(:,strcmp(header,'hip_adduction_r')),'-','Color',Cs); hold on;
    title('hip adduction');
    subplot(nr,nc,3);
    plot(R.Qs(:,strcmp(header,'hip_rotation_r')),'-','Color',Cs); hold on;
    title('hip rotation');
    
    
    subplot(nr,nc,5);
    plot(R.Qs(:,strcmp(header,'knee_angle_r')),'-','Color',Cs); hold on;
    title('knee angle');
    
    subplot(nr,nc,7);
    plot(R.Qs(:,strcmp(header,'ankle_angle_r')),'-','Color',Cs); hold on;
    title('ankle angle');
    subplot(nr,nc,8);
    plot(R.Qs(:,strcmp(header,'subtalar_angle_r')),'-','Color',Cs); hold on;
    title('subtalar angle');
    subplot(nr,nc,9);
    plot(R.Qs(:,strcmp(header,'mtp_angle_r')),'-','Color',Cs); hold on;
    title('mtp angle');
    
    for i=1:nc
        subplot(nr,nc,nr*nc-i+1)
        xlabel('% gait cycle');
    end
    for i=1:nr
        subplot(nr,nc,(i-1)*nc+1)
        ylabel('angle [deg]');
    end
    
    
    %% Plot other kinetics
    
    axes('parent', tab7);
    
    nr = 3; nc = 3;
    subplot(nr,nc,1);
    plot(R.Tid(:,strcmp(header,'hip_flexion_r')),'-','Color',Cs); hold on;
    title('hip flexion');
    subplot(nr,nc,2);
    plot(R.Tid(:,strcmp(header,'hip_adduction_r')),'-','Color',Cs); hold on;
    title('hip adduction');
    subplot(nr,nc,3);
    plot(R.Tid(:,strcmp(header,'hip_rotation_r')),'-','Color',Cs); hold on;
    title('hip rotation');
    
    
    subplot(nr,nc,5);
    plot(R.Tid(:,strcmp(header,'knee_angle_r')),'-','Color',Cs); hold on;
    title('knee angle');
    
    subplot(nr,nc,7);
    plot(R.Tid(:,strcmp(header,'ankle_angle_r')),'-','Color',Cs); hold on;
    title('ankle angle');
    subplot(nr,nc,8);
    plot(R.Tid(:,strcmp(header,'subtalar_angle_r')),'-','Color',Cs); hold on;
    title('subtalar angle');
    subplot(nr,nc,9);
    plot(R.Tid(:,strcmp(header,'mtp_angle_r')),'-','Color',Cs); hold on;
    title('mtp angle');
    
    
    for i=1:nc
        subplot(nr,nc,nr*nc-i+1)
        xlabel('% gait cycle');
    end
    for i=1:nr
        subplot(nr,nc,(i-1)*nc+1)
        ylabel('Joint moment [Nm]');
    end
    
    %% Ground reaction force  
     axes('parent', tab8);
    for i=1:6
        subplot(2,3,i);
        plot(R.GRFs(:,i),'-','Color',Cs); hold on;
        title(R.colheaders.GRF{i});
        xlabel('% stride');
    end
else
    warning(['File not found: ' ResultsFile]);
end

end

