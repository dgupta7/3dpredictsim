close all
clear
clc

% 
% bounds_mtj = [-35*0.75;35]*pi/180;
% 
% q_mt = linspace(-35,35,500)'*pi/180;
% 
% k_pass.mtj = [-30 5 5 -3]';
% theta.pass.mtj = bounds_mtj;
% tau_pass = k_pass.mtj(1,1)*exp(k_pass.mtj(2,1)*(q_mt-theta.pass.mtj(2,1))) + ...
%     k_pass.mtj(3,1)*exp(k_pass.mtj(4,1)*(q_mt-theta.pass.mtj(1,1)));
% 
% q0 = 5*pi/180;
% kMT_li = 10;
% M_li = -(q_mt - q0)*kMT_li;
% M_li = M_li + tau_pass;
% 
% figure
% plot(q_mt,M_li)
% grid on

%%
q_mt = linspace(0,20,500)'*pi/180;
% q_mt = 0;

%% long plantar ligament
calcnLPL2mtj = 0.053283;
mtj2mfLPL = 0.025085;
phi0_LPL = 0.96991;

a = calcnLPL2mtj;
b = mtj2mfLPL;
L0 = 0.0442/1.01;

phi = phi0_LPL + q_mt;
l = sqrt(a^2 + b^2 - 2*a*b.*cos(phi));
h = a*b./l.*sin(phi);

nu = 0.4;
A0 = 14.7*1.7;

lambda = l./L0;
A = A0*lambda.^(-nu*2);

% Gefen2001
a1 = -412640.5;
a2 = 2235967.7;
a3 = -4841544.8;
a4 = 5236972.7;
a5 = -2829945.7;
a6 = 611190.6;
s_LPL = a1*lambda.^5 + a2*lambda.^4 + a3*lambda.^3 + a4*lambda.^2 + a5*lambda + a6;
    
F_LPL = s_LPL.*A;

e = (lambda-1)*100;

M_LPL = -F_LPL.*h;

figure
plot(e,F_LPL)

%% short plantar ligament
calcnSPL2mtj = 0.025667;
mtj2mfSPL = 0.016521;
phi0_SPL = 1.0104;

a = calcnSPL2mtj;
b = mtj2mfSPL;
L0 = 0.0219/1.01;

phi = phi0_SPL + q_mt;
l = sqrt(a^2 + b^2 - 2*a*b*cos(phi));
h = a*b./l.*sin(phi);

nu = 0.4;
A0 = 11.6*1.4;

lambda = l./L0;
A = A0*lambda.^(-nu*2);

s_SPL = a1*lambda.^5 + a2*lambda.^4 + a3*lambda.^3 + a4*lambda.^2 + a5*lambda + a6;
F_SPL = s_SPL.*A;

e = (lambda-1)*100;

M_SPL = -F_SPL.*h;

figure
plot(e,F_SPL)

%%
figure
plot(q_mt*180/pi,M_LPL+M_SPL)
grid on
hold on


%%

x = linspace(-10,20,500)'*pi/180;

y = -8*(exp(4.5*(x+0.01))-1);
y2 = 8*(exp(-15*(x+0.15)));

plot(x*180/pi,y)
plot(x*180/pi,y2)
% plot(x*180/pi,y+y2)

%%

calcnPF2mtj = 0.06695;
mtj2mttPF = 0.091714;
phi0 = 2.1274;

q_mt = 0;
phi = phi0 + q_mt; % top angle of WL triangle
l_PF_fa = sqrt(calcnPF2mtj^2 + mtj2mttPF^2 - 2*calcnPF2mtj*mtj2mttPF*cos(phi)); % length of PF spanning arch
MA_PF = calcnPF2mtj*mtj2mttPF/l_PF_fa*sin(phi); % moment arm of PF to mtj


