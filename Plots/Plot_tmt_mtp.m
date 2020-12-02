function [] = Plot_tmt_mtp(ResultsFile)
if exist(ResultsFile,'file')
    load(ResultsFile,'R');
    
    [~,filename,~]= fileparts(ResultsFile);
    
    idx_tmt = strcmp(R.colheaders.joints  ,'tmt_angle_r');
    idx_mtp = strcmp(R.colheaders.joints  ,'mtp_angle_r');
    
    
    x = 1:(100-1)/(size(R.Qs,1)-1):100;
    x1 = R.Event.Stance;
    
    set(0,'defaultTextInterpreter','none');
    
    figure
    plot(x,R.Qs(:,idx_mtp),'DisplayName','mtp')
    hold on
    plot(x,R.Qs(:,idx_tmt)*10,'DisplayName','tmt (x10)') 
    line([x1 x1], get(gca, 'ylim'),'color','k','LineStyle','--','DisplayName','push-off');
    grid on
    ylabel('angle (°)')
    xlabel('Gait cycle (%)')
    title({'kinematics'; filename})
    legend('location','best')
    
    
    figure
    plot(x,R.Tid(:,idx_mtp),'DisplayName','mtp')
    hold on
    plot(x,R.Tid(:,idx_tmt),'DisplayName','tmt')
    line([x1 x1], get(gca, 'ylim'),'color','k','LineStyle','--','DisplayName','push-off');
    grid on
    ylabel('torque (Nm)')
    xlabel('Gait cycle (%)')
    title({'kinetics'; filename})
    legend('location','best')
    
    
else
    warning(['File not found: ' ResultsFile]);
end
end