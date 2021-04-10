% close all
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
q_mt = linspace(0,15,500)'*pi/180;
% q_mt = 0;

%% long plantar ligament
calcnLPL2mtj = 0.053283;
mtj2mfLPL = 0.025085;
phi0_LPL = 0.96991;

a = calcnLPL2mtj;
b = mtj2mfLPL;
L0 = 0.0442*1.02;

phi = phi0_LPL + q_mt;
l = sqrt(a^2 + b^2 - 2*a*b.*cos(phi));
h = a*b./l.*sin(phi);

nu = 0.4;
E = 260;
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

e = (lambda-1);

M_LPL = -F_LPL.*h;
M_L = -E*e.*A.*h;

% figure
% plot(e,F_LPL)

%% short plantar ligament
calcnSPL2mtj = 0.025667;
mtj2mfSPL = 0.016521;
phi0_SPL = 1.0104;

a = calcnSPL2mtj;
b = mtj2mfSPL;
L0 = 0.0219*1.02;

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
M_S = -E*e.*A.*h;

% figure
% plot(e,F_SPL)

%%
f1 = figure;
plot(q_mt*180/pi,M_LPL+M_SPL,'DisplayName','Gefen2001')
grid on
hold on
% plot(q_mt*180/pi,M_L+M_S,'DisplayName','E = 260MPa')
legend('Location','southwest')
xlabel('Angle (°)')
ylabel('Torque (Nm)')
title('Equivalent ligament stiffness midtarsal joint')
ylim([-100,100])

%%

x = linspace(-25,15,500)'*pi/180;

c1 = -8;
c2 = 8;
c3 = 0.2;
% y = c1*exp(c2*(x-c3));

y = -9*(exp(4*(x-2*pi/180))-1)*1.2;

y1 = 2*exp(-10*(x+0.1));

% plot(x*180/pi,y,'--')
plot(x*180/pi,y1)
y2 = (y+y1);
plot(x*180/pi,y2)

ylim([-20,20])

% figure
% plot(x,y+y2)
% xlim([-0.15,0.1])

%%

calcnPF2mtj = 0.06695;
mtj2mttPF = 0.091714;
phi0 = 2.1274;

q_mt = 0;
phi = phi0 + q_mt; % top angle of WL triangle
l_PF_fa = sqrt(calcnPF2mtj^2 + mtj2mttPF^2 - 2*calcnPF2mtj*mtj2mttPF*cos(phi)); % length of PF spanning arch
MA_PF = calcnPF2mtj*mtj2mttPF/l_PF_fa*sin(phi); % moment arm of PF to mtj


%%
k1 = 90;
dl_0 = 2*pi/180;

y = -k1*((x-0*pi/180) - dl_0*tanh((x-0*pi/180)/dl_0/1.25));
y2 = (-exp(20*(x-20*pi/180)) + exp(-25*(x+10*pi/180)))/2;

% plot(x*180/pi,y)
% plot(x*180/pi,y2)
% plot(x*180/pi,y+y2)


% interp1(x,(y+y2),0)


%%
k1 = 100;

y = -k1*x;
y2 = (exp(5*(x-10*pi/180)) + exp(-3*(x+20*pi/180)))/2;
% y2 = (exp(10*(x-15*pi/180)) + exp(-3*(x+20*pi/180)))/2;
% y2 = (exp(10*(x-15*pi/180)) + exp(-3*(x+20*pi/180)))/2;

% figure
% hold on
% plot(x*180/pi,y2)
% 
% 
% figure(f1)
% hold on
% grid on
% plot(x*180/pi,y)
% 
% plot(x*180/pi,y.*y2)
% ylim([-30,30])



%%

k1 = 50;
dl_0 = 1.5*pi/180;

% y = -k1*((x-1*pi/180) - dl_0*tanh((x-1*pi/180)/dl_0/1));

y = -1*exp(10*(x-0*pi/180)) + exp(-10*(x+0*pi/180));

% plot(x*180/pi,y)
% ylim([-30,30])


%%
kMT_li = 90;
q_mt = x;

M_li = -kMT_li*q_mt.*(exp(5*(q_mt-10*pi/180)) + exp(-3*(q_mt+20*pi/180)))/2;
y0 = M_li;


figure(f1)
hold on
grid on
% plot(x*180/pi,y0,'--','DisplayName','approx 1')
% plot(x*180/pi,-90*x,':','DisplayName','k = 90Nm/rad')

%%
M_li = (-exp(12*(q_mt-8*pi/180)) + 1.5*exp(-5*(q_mt+10*pi/180)))*5*3-2;

y = M_li;

% plot(x*180/pi,y,'DisplayName','approx 2')

interp1(y,x,0)*180/pi
interp1(x,y,0)


%%
sf_li = 0.2;
k1 = 12;
t1 = 8*pi/180*sf_li;
f2 = 1.5*(2-sf_li);
k2 = 5;
t2 = 10*pi/180*sf_li;
M_li = (-exp(k1*(q_mt-t1)) + f2*exp(-k2*(q_mt+t2)))*kMT_li*5/90;
y = M_li;
% plot(x*180/pi,y,'DisplayName',['approx 2 ' num2str(sf_li)])


%%

sf = 0.1;
t = 0;
c1 = 12;%*(2-t)/(1-t+sf);
c2 = 8*sf;

t = -0;
c3 = 1.5*(1+sf)/2;
c4 = 5;%*(2-t)/(1-t+sf);
c5 = 10*sf;

M_li = (-exp(c1*(q_mt-c2*pi/180)) + c3*exp(-c4*(q_mt+c5*pi/180)))*kMT_li*5/90;%*(1+sf)/2;

y = M_li;

% plot(x*180/pi,y,'DisplayName',['approx stiff ' num2str(sf)])


interp1(y,x,0)*180/pi
interp1(x,y,0)

%%

% c0 = 80;
% 
% c1 = 8;
% c2 = 10;
% 
% c3 = 120;
% 
% c4 = 4;
% c5 = 20;
% 
% M_li = -c0*exp(c1*(q_mt-c2*pi/180)) + c3*exp(-c4*(q_mt+c5*pi/180));

M_li = (-exp(12*(q_mt-3*pi/180)) + exp(-12*(q_mt+3*pi/180)))*10;


y = M_li;

% plot(x*180/pi,y,'DisplayName',['approx stiff 2'])

interp1(x,y,0)


%%
% ylim([-60,60])
%%

% xlim([-7,5])
% ylim([-5,5])



