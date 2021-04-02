

f_PF_stiffness = f_getPlantarFasciaStiffnessModelCasADiFunction('Gefen2001');


q_mtp = linspace(-45,45,500)*pi/180;

for i=1:length(q_mtp)
    [~,T_mtp(i)] = getPassiveMtjMomentWindlass_v3(0,0,q_mtp(i),f_PF_stiffness);
end

figure
plot(q_mtp*180/pi,T_mtp)
hold on

tau_mtp = 1*exp(-2*(q_mtp+0.2));

plot(q_mtp*180/pi,tau_mtp)
grid on
T = T_mtp + tau_mtp;
plot(q_mtp*180/pi,T)

% plot(q_mtp*180/pi,-17*q_mtp)