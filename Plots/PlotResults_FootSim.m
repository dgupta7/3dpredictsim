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
    plot(R.l_fa(i,:)*1000,R.Fs_tib,'Color',CsV(i,:))
    plot(R.l_fa_ext(i,:)*1000,R.Fs_tib,'--','Color',CsV(i,:))
    ylabel('tibia force (N)')
    xlabel('arch length(mm)')

    subplot(3,3,2)
    hold on
%     arch_compr = 1-(R.h_fa(i,:) - R.h0_fa(i,:))./R.h0_fa(i,:);
    plot(R.h_fa(i,:)*1000,R.Fs_tib,'Color',CsV(i,:))
    plot(R.h_fa_ext(i,:)*1000,R.Fs_tib,'--','Color',CsV(i,:))
    ylabel('tibia force (N)')
    xlabel('arch height (mm)')

%     subplot(3,3,3)
%     hold on
%     plot(R.Fs_tib,R.Qs(i,:,R.jointfi.tmt.r)*180/pi,'Color',CsV(i,:))
%     xlabel('tibia force (N)')
%     ylabel('tmt angle (°)')
    
    
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

%     figure(2)
%     xy = [squeeze(R.toes_or(i,1,1:2)),squeeze(R.metatarsi_or(i,1,1:2)),...
%         squeeze(R.calcn_or(i,1,1:2)),squeeze(R.talus_or(i,1,1:2)),...
%         squeeze(R.tibia_or(i,1,1:2))];
%     
%     hold on
%     grid on
%     plot(xy(1,:),xy(2,:),'-o','Color',CsV(i,:))
%     plot(xy(1,1),xy(2,1),'x','Color',CsV(i,:))
%     axis equal
    
end


CsV = hsv(n_tib);

leg1 = [];

figure
hold on
grid on
for i=1:n_tib
    
    
    xy = [squeeze(R.toes_or(1,i,1:2)),squeeze(R.metatarsi_or(1,i,1:2)),...
        squeeze(R.calcn_or(1,i,1:2)),squeeze(R.talus_or(1,i,1:2)),...
        squeeze(R.tibia_or(1,i,1:2))];
    
    z = [squeeze(R.toes_or(1,i,3)),squeeze(R.metatarsi_or(1,i,3)),...
        squeeze(R.calcn_or(1,i,3)),squeeze(R.talus_or(1,i,3)),...
        squeeze(R.tibia_or(1,i,3))];
    
    subplot(2,2,1)
    hold on
    plot(xy(1,4:5),xy(2,4:5),'-','Color',CsV(i,:))
    plot(xy(1,1),xy(2,1),'x','Color',CsV(i,:))
    plot(xy(1,2),xy(2,2),'o','Color',CsV(i,:))
    plot(xy(1,3),xy(2,3),'.','Color',CsV(i,:))
    plot(xy(1,4),xy(2,4),'*','Color',CsV(i,:))
    axis equal
    title('sagittal plane (right)')
    xlabel('x')
    ylabel('y')
    xlim([-0.07,0.12])
    ylim([0,0.1])
    
    subplot(2,2,2)
    hold on
    plot(-z(1,4:5),xy(2,4:5),'-','Color',CsV(i,:))
    p1=plot(-z(1,1),xy(2,1),'x','Color',CsV(i,:),'DisplayName',['F_y = ' num2str(R.Fs_tib(i)) 'N']);
    plot(-z(1,2),xy(2,2),'o','Color',CsV(i,:))
    plot(-z(1,3),xy(2,3),'.','Color',CsV(i,:))
    plot(-z(1,4),xy(2,4),'*','Color',CsV(i,:))
    leg1(i) = p1;
    axis equal
    title('frontal plane (front)')
    xlabel('-z')
    ylabel('y')
    xlim([-0.07,0.05])
    ylim([0,0.1])
    lg1=legend(leg1,'Location','northeast');
    subplot(2,2,3)
    hold on
    p6=plot(xy(1,4:5),-z(1,4:5),'-','Color',CsV(i,:),'DisplayName','tibia');
    p2=plot(xy(1,1),-z(1,1),'x','Color',CsV(i,:),'DisplayName','toes or');
    p3=plot(xy(1,2),-z(1,2),'o','Color',CsV(i,:),'DisplayName','metatarsi or');
    p4=plot(xy(1,3),-z(1,3),'.','Color',CsV(i,:),'DisplayName','calcn or');
    p5=plot(xy(1,4),-z(1,4),'*','Color',CsV(i,:),'DisplayName','talus or');
    leg2 = [p2,p3,p4,p5,p6];
    lg2=legend(leg2,'Location','northeast');
    axis equal
    title('transverse plane (top)')
    xlabel('x')
    ylabel('-z')
    xlim([-0.07,0.12])
    ylim([-0.07,0.07])
    
end
title(lg1,'colour meaning')
lhPos = lg1.Position;
lhPos(2) = lhPos(2)-0.4;
lhPos(1) = lhPos(1)+0.1;
set(lg1,'position',lhPos);

title(lg2,'symbol meaning')
lhPos = lg2.Position;
lhPos(1) = lhPos(1)+0.3;
set(lg2,'position',lhPos);

end