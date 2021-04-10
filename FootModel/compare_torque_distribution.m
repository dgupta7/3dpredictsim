clear
clc
AddCasadiPaths();

S.sf_PF = 1;
S.sf_li = 1;
PF_stiffness = {'Gefen2001'};
% PF_stiffness = {'Natali2010'};
% PF_stiffness = {'linear'};

ls = 0.135;
f_PF_stiffness = f_getPlantarFasciaStiffnessModelCasADiFunction(PF_stiffness,'ls',ls);

S.dMT = 0;
S.MT_li_nonl = 1;
S.kMT_li = 90;

q_mt = linspace(-20,20,1000)'*pi/180;
q_mtp = [-30,0,30]*pi/180;
% q_mtp = 0;

Fmin = 1;
Fmax = 30;

%%

[M0,~,M1,M2,~,~,~,~,~,~,~] = ...
                        getPassiveMtjMomentWindlass_v3(5*pi/180,0,0*pi/180,f_PF_stiffness,S);

%%
for j=1:length(q_mtp)
    for i=1:length(q_mt)
        [M0,~,M1,M2,~,~,~,~,~,~,~] = ...
                        getPassiveMtjMomentWindlass_v3(q_mt(i),0,q_mtp(j),f_PF_stiffness,S);
        M_PF(i,j) = M1;
        M_li(i,j) = M2;
        M(i,j) = M0;
        
        M_PF_rel(i,j) = M1/M2;
        
        if abs(M_PF_rel(i,j))<1e-2
            M_PF_rel(i,j) = 0;
        end
    end
    k(:,1) = -(M(2:end,j)-M(1:end-1,j))./(q_mt(2:end)-q_mt(1:end-1));
    k_eq(:,j) = [k(1);k];
    
    
    idx = find(M(:,j)<=-Fmin & M(:,j)>=-Fmax);
    
%     E = trapz(q_mt(idx),M(idx));
%     Ep_tot(1,j) = E;
%     
%     E = trapz(q_mt(idx),M_PF(idx));
%     Ep_PF(1,j) = E;
%     
%     E = trapz(q_mt(idx),M_li(idx));
%     Ep_li(1,j) = E;
    
    temp1 = zeros(length(idx),1);
    temp2 = zeros(length(idx),1);
    temp3 = zeros(length(idx),1);
    for ii=2:length(idx)
       temp1(ii) = trapz(q_mt(1:idx(ii)),M(1:idx(ii),j));
       temp2(ii) = trapz(q_mt(1:idx(ii)),M_PF(1:idx(ii),j));
       temp3(ii) = trapz(q_mt(1:idx(ii)),M_li(1:idx(ii),j));
    end
    temp1(2:end) = temp1(2:end)-temp1(2);
    temp2(2:end) = temp2(2:end)-temp2(2);
    temp3(2:end) = temp3(2:end)-temp3(2);
    
    Ep_f.(['q' num2str(j)]).E_t = temp1;
    Ep_f.(['q' num2str(j)]).E_PF = temp2;
    Ep_f.(['q' num2str(j)]).E_li = temp3;
    Ep_f.(['q' num2str(j)]).q_mt = q_mt(idx);
    Ep_f.(['q' num2str(j)]).M = M(idx,j);
    Ep_tot(1,j) = temp1(end);
    Ep_PF(1,j) = temp2(end);
    Ep_li(1,j) = temp3(end);
    
    q_t0(1,j) = interp1(M(:,j),q_mt,0);
    q_t0(2,j) = interp1(q_mt,M_li(:,j),q_t0(1,j));
    
end


figure
subplot(2,3,1)
for j=1:length(q_mtp)
    hold on
    grid on
    plot(q_mt*180/pi,M(:,j),'DisplayName',['q_{mtp} = ' num2str(q_mtp(j)*180/pi)])
    
end
xlabel('mtj angle (°)')
ylabel('mtj torque (Nm)')
legend('Location','best');
title('Total midtarsal joint torque')
ylim([-100,30])

subplot(2,3,2)
for j=1:length(q_mtp)
    hold on
    grid on
    p1=plot(q_mt*180/pi,M_PF(:,j),'DisplayName',['PF, q_{mtp} = ' num2str(q_mtp(j)*180/pi)]);
    ps(j)=p1;
    
end
p2=plot(q_mt*180/pi,M_li(:,1),'DisplayName','li');
ps(j+1)=p2;
for j=1:length(q_mtp)
    hold on
    grid on
    plot(q_t0(1,j)*180/pi,q_t0(2,j),'v','Color',ps(j).Color)
    
end
xlabel('mtj angle (°)')
ylabel('mtj torque (Nm)')
legend(ps,'Location','best');
title('Midtarsal joint torque components')
ylim([-100,30])

subplot(2,3,3)
for j=1:length(q_mtp)
    plot(q_mt*180/pi,db(M_PF_rel(:,j)),'DisplayName',['q_{mtp} = ' num2str(q_mtp(j)*180/pi)])
    hold on
    grid on
    
end
xlabel('mtj angle (°)')
ylabel('PF/PL torque (dB)')
legend('Location','best');
title('Relative PF torque wrt plantar ligaments')

subplot(2,3,4)
for j=1:length(q_mtp)
    plot(q_mt*180/pi,k_eq(:,j),'DisplayName',['q_{mtp} = ' num2str(q_mtp(j)*180/pi)])
    hold on
    grid on
    
end
xlabel('mtj angle (°)')
ylabel('k (Nm/rad)')
legend('Location','best');
title('Midtarsal stiffness')
    
subplot(2,3,5)
for j=1:length(q_mtp)
    h1=plot(-Ep_f.(['q' num2str(j)]).M,-Ep_f.(['q' num2str(j)]).E_t,...
        'DisplayName',['q_{mtp} = ' num2str(q_mtp(j)*180/pi)]);
    hold on
    grid on
    hs(j)=h1;
    
end
% for j=1:length(q_mtp)
%     plot(-Ep_f.(['q' num2str(j)]).M,-Ep_f.(['q' num2str(j)]).E_PF,...
%         '--','Color',hs(j).Color','DisplayName','PF');
%     plot(-Ep_f.(['q' num2str(j)]).M,-Ep_f.(['q' num2str(j)]).E_li,...
%         ':','Color',hs(j).Color','DisplayName','LPL + SPL');
% end
xlabel('External torque (Nm)')
ylabel('Potential energy (J)')
legend('Location','best');
title('Cumulative energy storage')

% subplot(2,3,5)
% for j=1:length(q_mtp)
%     h1=plot(Ep_f.(['q' num2str(j)]).q_mt*180/pi,-Ep_f.(['q' num2str(j)]).E_t,...
%         'DisplayName',['q_{mtp} = ' num2str(q_mtp(j)*180/pi)]);
%     hold on
%     grid on
% 
% end
% xlabel({'Midtarsal angle (°)',['for T_{ext} = ' num2str(Fmin) '..' num2str(Fmax) 'Nm']})
% ylabel('Potential energy (J)')
% legend('Location','best');
% title('Cumulative energy storage')


subplot(2,3,6)
% plot(q_mtp*180/pi,-Ep_tot,'o','DisplayName','Total')
% hold on
% plot(q_mtp*180/pi,-Ep_PF,'d','DisplayName','PF')
% plot(q_mtp*180/pi,-Ep_li,'^','DisplayName','LPL + SPL')

Ep = [-Ep_tot',-Ep_PF',-Ep_li']';

for j=1:length(q_mtp)
    leg{j} = ['q_{mtp} = ' num2str(q_mtp(j)*180/pi) '°'];
end

c = categorical({'Total','Plantar fascia','Other ligaments'});
c = reordercats(c,[3;2;1]);

bar(c,Ep)
ylabel('Potential energy (J)')
legend(leg,'Location','northeast');
title('Total energy storage')










