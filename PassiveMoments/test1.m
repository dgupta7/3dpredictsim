
itmt = find(strcmp(R.colheaders.joints,'tmt_angle_r'));
imtp = find(strcmp(R.colheaders.joints,'mtp_angle_r'));
q_tmt = R.Qs(:,itmt);
qdot_tmt = R.Qdots(:,itmt);
q_mtp = R.Qs(:,imtp);


kTMT_li = 1.5/(pi/180)/5;
kTMT_PF = R.S.kTMT;
dTMT = R.S.dTMT;
cWL = 0.03;

x = 1:(100-1)/(size(R.Qs,1)-1):100;

for i=1:length(R.Qs)
    [Mi, M_PFi,F_PFi,~,~,li,l0i,L0,hi,h0i,H0] = ...
        getPassiveTmtjMomentWindlass(q_tmt(i)*pi/180,qdot_tmt(i),q_mtp(i)*pi/180,kTMT_li,kTMT_PF,dTMT,R.S.subject,cWL);
    
    M(i) = Mi;
    M_PF(i) = M_PFi;
    l(i) = li;
    h(i) = hi;
    F_PF(i) = F_PFi;
    l0(i) = l0i;
    h0(i) = h0i;
end


label_fontsize  = 12;

figure
subplot(2,3,1)
hold on
plot(x,q_tmt)
title('tmt angle')
xlabel('Gait cycle (%)','Fontsize',label_fontsize);
ylabel('Angle (°)','Fontsize',label_fontsize);

subplot(2,3,4)
hold on
plot(x,q_mtp)
title('mtp angle')
xlabel('Gait cycle (%)','Fontsize',label_fontsize);
ylabel('Angle (°)','Fontsize',label_fontsize);

subplot(2,3,2)
hold on
plot(x,R.Tid(:,itmt),'DisplayName','T')
plot(x,M,'DisplayName','M')
plot(x,-M_PF,'DisplayName','M PF')
legend('location','best')
title('tmt moment')
xlabel('Gait cycle (%)','Fontsize',label_fontsize);
ylabel('Moment (Nm)','Fontsize',label_fontsize);

subplot(2,3,5)
hold on
plot(x,F_PF)
title('Plantar fascia force')
xlabel('Gait cycle (%)','Fontsize',label_fontsize);
ylabel('Force (N)','Fontsize',label_fontsize);

subplot(2,3,3)
hold on
plot(x,L0*ones(size(x)),'DisplayName','unloaded')
plot(x,l0,'DisplayName','Windlass')
plot(x,l,'DisplayName','PF stretched')
legend('location','best')
title('Foot arch length')
xlabel('Gait cycle (%)','Fontsize',label_fontsize);
ylabel('relative length (-)','Fontsize',label_fontsize);

subplot(2,3,6)
hold on
plot(x,H0*ones(size(x)),'DisplayName','unloaded')
plot(x,h0,'DisplayName','Windlass')
plot(x,h,'DisplayName','PF stretched')
legend('location','best')
title('Foot arch height')
xlabel('Gait cycle (%)','Fontsize',label_fontsize);
ylabel('relative height (-)','Fontsize',label_fontsize);



