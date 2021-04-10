


f_PF_stiffness = f_getPlantarFasciaStiffnessModelCasADiFunction('Barrett2018');
% f_PF_stiffness = f_getPlantarFasciaStiffnessModelCasADiFunction('Natali2010');


q_mtp = linspace(-45,45,500)*pi/180;

for i=1:length(q_mtp)
    [~,T_mtp(i)] = getPassiveMtjMomentWindlass_v3(0,0,q_mtp(i),f_PF_stiffness);
end

% old model
T_pass_mtp = -1.5/(pi/180)/5.*q_mtp;


% plot(q_mtp*180/pi,T_mtp)

% % k1 = 1;
% % dq = 0.2;
% k1 = 5;
% dq = 0.2;
% tau_mtp = k1*exp(-2*(q_mtp+dq));

% plot(q_mtp*180/pi,tau_mtp)
% grid on
% T = T_mtp + tau_mtp;
% plot(q_mtp*180/pi,T)

figure
hold on
grid on
plot(q_mtp,T_pass_mtp-T_mtp)

%%
% T_comp = 4 - 10.*q_mtp; % linear
% T_comp = 9*q_mtp.^2 - 5*q_mtp + 5; % Natali
% T_comp = 1 - 15.*q_mtp; % Gefen
% T_comp = 4 - 12.*q_mtp; % Cheng
T_comp =  -16.*q_mtp; % Barrett


plot(q_mtp,T_comp)

%%

figure
plot(q_mtp,T_comp+T_mtp)
hold on
grid on
plot(q_mtp,T_pass_mtp)

