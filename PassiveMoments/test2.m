
q_tmt = (-15:1:15)*pi/180;
q_mtp = (-45:15:45)*pi/180;
qdot_tmt = 0;
kTMT_li = 17;
kTMT_PF = 800;
dTMT = 0;
subject = 's1_Poggensee';
cWL = 0.03;

M = zeros(length(q_tmt),length(q_mtp));

figure
hold on
grid on
for j=1:length(q_mtp)
    for i=1:length(q_tmt)
   
        M(i,j) = getPassiveTmtjMomentWindlass(q_tmt(i),qdot_tmt,q_mtp(j),kTMT_li,kTMT_PF,dTMT,subject,cWL);

    end
    plot(q_tmt*180/pi,M(:,j),'DisplayName',num2str(q_mtp(j)*180/pi))
end
legend
