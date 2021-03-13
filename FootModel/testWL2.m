
% AddCasadiPaths();
% import casadi.*
% % f_AllPassiveTorques = Function.load('f_AllPassiveTorques');
% f_passiveWLTorques = Function.load('f_passiveWLTorques');




q_tmt = (-15:1:15)*pi/180;
q_mtp = [-45,0]*pi/180;
qdot_tmt = 0;
kTMT_li = 17;
kTMT_PF = 2000;
dTMT = 0;
subject = 's1_Poggensee';
cWL = 0.03;

M = zeros(length(q_tmt),length(q_mtp));
M1 = zeros(length(q_tmt),length(q_mtp));

figure
hold on
grid on
for j=1:length(q_mtp)
    for i=1:length(q_tmt)
   
        M(i,j) = getPassiveTmtjMomentWindlass(q_tmt(i),qdot_tmt,q_mtp(j),kTMT_li,kTMT_PF,dTMT,subject,cWL);
        M(i,j) = getPassiveTmtjMomentWindlass1(q_tmt(i),qdot_tmt,q_mtp(j),kTMT_li,kTMT_PF,dTMT,subject,cWL);
%         M1(i,j) = full(f_passiveWLTorques(q_tmt(i),qdot_tmt,q_mtp(j)));

    end
    plot(q_tmt*180/pi,M(:,j),'DisplayName',num2str(q_mtp(j)*180/pi))
    plot(q_tmt*180/pi,M1(:,j),'--o','DisplayName',num2str(q_mtp(j)*180/pi))
end
leg = legend('location','best');
title(leg,'q mtp')
xlabel('q tmt')
ylabel('T tmt')
