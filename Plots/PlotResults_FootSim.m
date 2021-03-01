function [] = PlotResults_FootSim(R)

n_mtp = length(R.Qs_mtp);
n_tib = length(R.Fs_tib);

CsV = hsv(n_mtp);

% h_fa_max = max(max(R.h_fa));
% h_fa_ext_max = max(max(R.h_fa_ext));

for i=1:n_mtp
    
    figure(1)
    
    subplot(3,3,1)
    hold on
    plot(R.l_fa(i,:)-R.L0,R.Fs_tib,'Color',CsV(i,:))
    plot(R.l_fa_ext(i,:)-R.L0,R.Fs_tib,'--','Color',CsV(i,:))
    ylabel('tibia force (N)')
    xlabel('arch length increase(mm)')

    subplot(3,3,2)
    hold on
%     arch_compr = 1-(R.h_fa(i,:) - R.h0_fa(i,:))./R.h0_fa(i,:);
    plot(R.h_fa(i,:),R.Fs_tib,'Color',CsV(i,:))
    plot(R.h_fa_ext(i,:),R.Fs_tib,'--','Color',CsV(i,:))
    ylabel('tibia force (N)')
    xlabel('arch height compression (-)')

    subplot(3,3,3)
    hold on
    plot(R.Fs_tib,R.Qs(i,:,R.jointfi.tmt.r)*180/pi,'Color',CsV(i,:))
    xlabel('tibia force (N)')
    ylabel('tmt angle (°)')
    
    
    subplot(3,3,4)
    hold on
    plot(R.Fs_tib,R.Qs(i,:,R.jointfi.tmt.r)*180/pi,'Color',CsV(i,:))
    xlabel('tibia force (N)')
    ylabel('tmt angle (°)')
    
    subplot(3,3,5)
    hold on
    plot(R.Fs_tib,R.Qs(i,:,R.jointfi.ankle.r)*180/pi,'Color',CsV(i,:))
    xlabel('tibia force (N)')
    ylabel('ankle angle (°)')

    subplot(3,3,6)
    hold on
    plot(R.Fs_tib,R.Qs(i,:,R.jointfi.subt.r)*180/pi,'Color',CsV(i,:),'DisplayName',['mtp: ' num2str(R.Qs_mtp(i)*180/pi) '°'])
    xlabel('tibia force (N)')
    ylabel('subt angle (°)')
    legend('Location','best')
    
    
    subplot(3,3,7)
    hold on
    plot(R.Fs_tib,R.GRF_calcn(i,:,2),'Color',CsV(i,:))
    xlabel('tibia force (N)')
    ylabel('GRF_y calcn (N)')

    subplot(3,3,8)
    hold on
    plot(R.Fs_tib,R.GRF_metatarsi(i,:,2),'Color',CsV(i,:))
    xlabel('tibia force (N)')
    ylabel('GRF_y metatarsi (N)')
    
    subplot(3,3,9)
    hold on
    plot(R.Fs_tib,R.F_PF(i,:),'Color',CsV(i,:))
    xlabel('tibia force (N)')
    ylabel('F PF (N)')

    figure(2)
    xy = [squeeze(R.toes_or(i,1,1:2)),squeeze(R.metatarsi_or(i,1,1:2)),...
        squeeze(R.calcn_or(i,1,1:2)),squeeze(R.talus_or(i,1,1:2)),...
        squeeze(R.tibia_or(i,1,1:2))];
    
    hold on
    grid on
    plot(xy(1,:),xy(2,:),'-o','Color',CsV(i,:))
    plot(xy(1,1),xy(2,1),'x','Color',CsV(i,:))
    axis equal
    
end


CsV = hsv(n_tib);

figure
hold on
grid on
for i=1:n_tib
    
    
    xy = [squeeze(R.toes_or(1,i,1:2)),squeeze(R.metatarsi_or(1,i,1:2)),...
        squeeze(R.calcn_or(1,i,1:2)),squeeze(R.talus_or(1,i,1:2)),...
        squeeze(R.tibia_or(1,i,1:2))];
    
    
    plot(xy(1,:),xy(2,:),'-o','Color',CsV(i,:))
    plot(xy(1,1),xy(2,1),'x','Color',CsV(i,:))
    axis equal
    
    
    
end


end