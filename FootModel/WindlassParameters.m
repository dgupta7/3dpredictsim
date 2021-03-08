%% Windlass mechanism
clear all
close all
clc

a = 0.1140;    % calcn length
b = 0.0666;    % metatarsi length
phi0 = 2.2114;    % tmt vector angle
g0 = phi0;      % tmt vector angle inf stiff

q1 = [-45:0.2:45]'*pi/180; % possible tmt angles
q2 = [-30:30:30]'*pi/180; % possible mtp angles

colr = hsv(length(q2));
L0 = sqrt(a^2 + b^2 - 2*a*b*cos(phi0));
H0 = 0.0374;

kTMT_l = 1000;
kTMT_li = 1.5/(pi/180)/5; %same as mtp, stiffness of ligaments and soft tissue not including PF

% weight factors to change cWL
% ws = [0.8:0.05:1.25];
ws = 1;
% ws = [0.005,0.01,0.02,0.04]/0.03;

for w=1:length(ws)

%% Full geometric model
% for infinitely stiff tendon
cWL = 0.03/0.97 * ws(w);

l0 = (1-cWL*(q2*180/pi)/20)*L0;
h0 = a*b./l0*sin(phi0);

g = acos( (a^2 + b^2 - l0.^2)/(2*a*b) );
q1_0 = g-g0;

% for elastic PF
phi = phi0 + q1;
l = sqrt(a^2 + b^2 - 2*a*b*cos(phi));
h = a*b./l.*sin(phi);


% get plantar fascia k
q1k = [-5:1:5]*pi/180;
M1 = kTMT_l*q1k;

lk = sqrt(a^2 + b^2 - 2*a*b*cos(phi0 + q1k));
hk = a*b./lk.*sin(phi0 + q1k);

k0s1 = M1./(hk.*(lk-L0));
k0s2 = M1./(H0*(lk-L0));
k01 = nanmean(k0s1);
k02 = nanmean(k0s2);

k0 = k01;
k2 = k0 * 120;

k0s_nl = kTMT_l*q1./(h.*(l-L0));
k0_nl = nanmean(k0s_nl);

% values from doi:10.1016/j.clinbiomech.2004.06.002
E_PF = 350; %MPa
A_PF = 290; %mm2
k_PF = E_PF*A_PF/L0;

% figure
% plot(q1k,k0s1)
% hold on
% grid on
% plot(q1k,k01*ones(size(q1k)))
% plot(q1k,k0s2)
% plot(q1k,k02*ones(size(q1k)))
% plot(q1k,k_PF*ones(size(q1k)))


% full tmt moment
M_li = kTMT_li*q1;
for i=1:length(q2)
    dl(:,i) = l(:)-l0(i);
    dl_pos(:,i) = dl(:,i) .*(tanh(dl(:,i)*1e6)+1)/2;
    F(:,i) = k0*dl_pos(:,i);
    F2(:,i) = k2*(dl_pos(:,i)).^2;
    M_te(:,i) = F(:,i).*h(:);
    M_te2(:,i) = F2(:,i).*h(:);
    M_f(:,i) = M_te(:,i) + M_li(:,1);
    M_f2(:,i) = M_te2(:,i) + M_li(:,1);
    c1(:,i) = dl_pos(:,i).*h(:);
    k_rot(:,i) = M_te(:,i)./(q1(:)-q1_0(i));
    k_rot2(:,i) = M_te2(:,i)./(q1(:)-q1_0(i));
%     k_rot(:,i) = M_te(:,i)./q1(:);
%     k_rot2(:,i) = M_te2(:,i)./q1(:);
end

% f1=figure;
% plot(dl(:,3),dl_pos(:,3))
% xlabel('\Delta l')
% ylabel('\Delta l > 0')
% grid on
% axis equal

% f2=figure;
% hold on
% grid on
% for j=1:i
%     plot(q1*180/pi,M_f(:,j),'color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi) '°'])
% end
% leg=legend('Location','best');
% xlabel('tmt angle (°)')
% ylabel('tmt moment (Nm)')
% title('Windlass mechanism')
% title(leg,'mtp angle')

% figure;
% hold on
% grid on
% for j=1:i
%     plot(q1*180/pi,M_f2(:,j),'color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi) '°'])
% end
% leg=legend('Location','best');
% xlabel('tmt angle (°)')
% ylabel('tmt moment (Nm)')
% title('Windlass mechanism (k2)')
% title(leg,'mtp angle')

% figure;
% hold on
% grid on
% for j=1:i
%     plot((q1(:)-q1_0(j))*180/pi,k_rot(:,j),'color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi) '°'])
%     plot((q1(:)-q1_0(j))*180/pi,k_rot2(:,j),'--','color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi) '°'])
% end
% leg=legend('Location','best');
% xlabel('tmt angle change from load(°)')
% ylabel('k (Nm/rad)')
% title('torsion stiffness')
% title(leg,'mtp angle')


% f3=figure;
% hold on
% grid on
% for j=1:i
%     plot((q1-q1_0(j))*180/pi,dl(:,j),'color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi) ' full'])
% end
% legend
% ylabel('\Delta l')
% xlabel('\Delta q1')
% title('Effect WL on l')


%% Approximations
% q1_0
cWLq = nanmean(q1_0./(q2*cWL));
% cWLq = -12.3342;
q1_0_lin = cWL*cWLq*q2;

% cWLq = nanmean(q1_0./(q2));
% q1_0_lin = cWLq*q2;

% figure
% plot(q2*180/pi,q1_0*180/pi)
% hold on
% grid on
% plot(q2*180/pi,q1_0_lin*180/pi,'--')
% xlabel('mtp angle (°)')
% ylabel('tmt angle inf stiff (°)')

% l
cWLl = nanmean((l/L0-1)./q1);
% cWLl = 0.2293;
l_lin = L0*(1 + cWLl*(q1+cos(q1)-1) );

% figure
% plot(q1*180/pi,l/L0)
% hold on
% grid on
% plot(q1*180/pi,l_lin/L0,'--')
% xlabel('tmt angle (°)')
% ylabel('normalized arch length (-)')

% h
cWLh = nanmean((h(q1~=0)-H0)./(q1(q1~=0)));
% cWLh = -0.0364;
h_lin = cWLh*q1+H0;

% figure
% plot(q1*180/pi,h/H0)
% hold on
% grid on
% plot(q1*180/pi,h_lin/H0,'--')
% xlabel('tmt angle (°)')
% ylabel('normalized arch height (-)')

% figure
% plot(q1*180/pi,h*1e3)
% hold on
% grid on
% plot(q1*180/pi,h_lin*1e3,'--')
% xlabel('tmt angle (°)')
% ylabel('h (mm)')
% title('arch height')

% dl
for i=1:length(q2)
    dl_lin(:,i) = l_lin(:)-l0(i);
end

% figure(f3)
% for j=1:i
%     plot((q1-q1_0_lin(j))*180/pi,dl_lin(:,j),'--','color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi) ' q1,l lin'])
% end


% lin tmt

for i=1:length(q2)
    dl_lin(:,i) = l_lin(:)-l0(i);
    dl_pos_lin(:,i) = dl_lin(:,i) .*(tanh(dl_lin(:,i)*1e6)+1)/2;
    F_lin(:,i) = k0*dl_pos_lin(:,i);
    M_te_lin(:,i) = F_lin(:,i).*h_lin(:);
    M_f_lin(:,i) = M_te_lin(:,i) + M_li(:,1);
    c1_lin(:,i) = dl_pos_lin(:,i).*h_lin(:);
end

% figure(f2)
% for j=1:i
%     plot(q1*180/pi,M_f_lin(:,j),'--','color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi) ' lin'])
% end

dc1 = (c1-c1_lin)./c1;
M_rel = (M_f-M_f_lin)./M_f;

% figure
% hold on
% grid on
% for j=1:i
%     plot(q1*180/pi,M_rel(:,j),'--','color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi)])
% end
% xlabel('tmt angle (°)')
% ylabel('relative linearization error tmt moment')


%
cw(1,w) = cWL;
cw(2,w) = cWLq;
cw(3,w) = cWLl;
cw(4,w) = cWLh;

end

%% influence of cWL on other linearization parameters
if w>1
    
    figure
    subplot(411)
    plot(ws,cw(1,:))
    title('cWL')
    
    subplot(412)
    plot(ws,cw(2,:))
    cm2 = nanmean(cw(2,:));
    hold on
    line(get(gca, 'xlim'),[1 1]*cm2,'color','k','LineStyle','--')
    title('cWLq')
    
    subplot(413)
    plot(ws,cw(3,:))
    cm3 = nanmean(cw(3,:));
    hold on
    line(get(gca, 'xlim'),[1 1]*cm3,'color','k','LineStyle','--')
    title('cWLl')
    
    subplot(414)
    plot(ws,cw(4,:))
    cm4 = nanmean(cw(4,:));
    hold on
    line(get(gca, 'xlim'),[1 1]*cm4,'color','k','LineStyle','--')
    title('cWLh')
    
    
    
end

%% more plots

% figure
% plot(q2*180/pi,h0*1e3)
% 
% figure
% plot(q2*180/pi,l0*1e3)

% figure
% hold on
% for j=1:i
%     plot(q1*180/pi,dl(:,j)*1000,'color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi)])
% end
% % plot(q1_0*180/pi,(l0-L0)*1000,'-k','DisplayName','l0-L0')
% line([-15,15],[0,0],'color','k')
% legend('Location','best')
% ylabel('horizontal elongation (mm)')
% xlabel('tmt angle (°)')

%% Linking to stiffness defined as vertical force / arch elongation

q_tmt = q1;
phi = phi0 + q_tmt;
l = sqrt(a^2 + b^2 - 2*a*b*cos(phi));
h = a*b./l.*sin(phi);


l1 = sqrt(a^2-h.^2);
l2 = sqrt(b^2-h.^2);

% l1+l2-l

k2 = k0 * 120;
BW = 64*9.81;

for i=1:length(q2)
    F2(:,i) = k2*(dl_pos(:,i)).^2;
    M_te2(:,i) = F2(:,i).*h(:);
    F_v1(:,i) = M_te(:,i).*(l1+l2)./(l1.*l2);
    F_v2(:,i) = M_te2(:,i).*(l1+l2)./(l1.*l2);
%     arch_compr1(:,i) = h(find(F_v1(:,i)<=0,1,'last')) - h(:);
%     arch_compr2(:,i) = h(find(F_v2(:,i)<=0,1,'last')) - h(:);
%     ac1_rel(:,i) = arch_compr1(:,i)./max(arch_compr1(:,i));
%     ac2_rel(:,i) = arch_compr2(:,i)./max(arch_compr2(:,i));
end

% ac1_rel = arch_compr1;
% ac2_rel = arch_compr2;

% % index 1 (q_mtp = -30°)
% idx = 1;
% idx_ac1 = find(F_v1(:,idx)<=BW);
% F_v1_ac1 = F_v1(idx_ac1,idx);
% ac1 = h0(idx) - h(idx_ac1);
% ac1_rel = ac1./max(ac1);



figure
hold on
grid on
for j=1:i
%     idx_ac1 = find(F_v1(:,j)<=BW);
%     F_v1_ac1 = F_v1(idx_ac1,j);
%     h0_ac1 = h(find(F_v1(:,j)<=0,1,'last'));
%     ac1 = h0_ac1 - h(idx_ac1);
%     ac1_max(j) = max(ac1);
%     ac1_rel = ac1./ac1_max(1);
%     plot(ac1_rel(:),F_v1_ac1(:)/BW,'--','color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi)])
    
    idx_ac2 = find(F_v2(:,j)<=BW);
    F_v2_ac2 = F_v2(idx_ac2,j);
    h0_ac2 = h(find(F_v2(:,j)<=0,1,'last'));
    ac2 = h0_ac2 - h(idx_ac2);
    ac2_max(j) = max(ac2);
    ac2_rel = ac2./ac2_max(1);
    plot(ac2_rel(:),F_v2_ac2(:)/BW,'-','color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi)])
    
end
xlabel('arch compression (-)')
ylabel('vertical force /BW (-)')
title({'arch stiffness','defined as by Welte et al, 2018'})
leg = legend('Location','best');
title(leg,'mtp angle')
xlim([0,1])


% figure
% hold on
% grid on
% for j=1:i
%     plot((l-L0)*1000,F_v1(:,j)/1000,'--','color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi)])
%     plot((l-L0)*1000,F_v2(:,j)/1000,'-','color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi)])
% end
% xlabel('horizontal elongation (mm)')
% ylabel('vertical force (kN)')
% title({'arch stiffness','defined as by Ker et al, 1987'})
% legend('Location','best')
% ylim([-0.5,4.5]);

% figure
% hold on
% grid on
% for j=1:i
%     plot((H0-h)*1000,F_v1(:,j)/1000,'--','color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi)])
%     plot((H0-h)*1000,F_v2(:,j)/1000,'-','color',colr(j,:),'DisplayName',[num2str(q2(j)*180/pi)])
% end
% xlabel('vertical displacement (mm)')
% ylabel('vertical force (kN)')
% title({'arch stiffness','defined as by Stearne et al, 2016'})
% legend('Location','best')
% ylim([-0.5,4.5]);

%%
% q_tmt = q1(q1>=0);
% phi = phi0 + q_tmt;
% l = sqrt(a^2 + b^2 - 2*a*b*cos(phi));
% h = a*b./l.*sin(phi);
% 
% % alpha = asin(a/h);
% % beta = asin(b/h);
% 
% l1 = sqrt(a^2-h.^2);
% l2 = sqrt(b^2-h.^2);
% 
% % l1+l2-l
% 
% k2 = k0 * 120;
% 
% F_vert = h.*(l-L0)*k0.*(l1+l2)./(l1.*l2);
% F_vert2 = h.*(l-L0).^2*k2.*(l1+l2)./(l1.*l2);

% figure
% plot((l-L0)*1000,F_vert/1000)
% hold on
% plot((l-L0)*1000,F_vert2/1000)
% xlabel('horizontal elongation (mm)')
% ylabel('vertical force (kN)')
% title({'arch stiffness','defined as by Ker et al, 1987'})
% legend('F = k \Deltal','F = k (\Deltal)^2')
% 
% figure
% plot((H0-h)*1000,F_vert/1000)
% hold on
% plot((H0-h)*1000,F_vert2/1000)
% xlabel('vertical displacement (mm)')
% ylabel('vertical force (kN)')
% title({'arch stiffness','defined as by Stearne et al, 2016'})
% legend('F = k \Deltal','F = k (\Deltal)^2')
% 
