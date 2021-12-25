

n = 500;

q_mtj = linspace(-20,30,n)*pi/180;



figure
hold on

%%
% 
% K_pass = [-50 20 60 -30]';
% theta_pass = [-0.4 0.3]';
% 
% Tau_pass_1 = K_pass(1,1)*exp(K_pass(2,1)*(q_mtj-theta_pass(2,1)));
% Tau_pass_2 = K_pass(3,1)*exp(K_pass(4,1)*(q_mtj-theta_pass(1,1)));
% 
% plot(q_mtj*180/pi,Tau_pass_1+Tau_pass_2,'--')

%%
K_pass = [-50 15 60 -30]';
theta_pass = [-0.4 0.5]';

Tau_pass_1 = K_pass(1,1)*exp(K_pass(2,1)*(q_mtj-theta_pass(2,1)));
Tau_pass_2 = K_pass(3,1)*exp(K_pass(4,1)*(q_mtj-theta_pass(1,1)));

plot(q_mtj*180/pi,Tau_pass_1+Tau_pass_2,'--')


