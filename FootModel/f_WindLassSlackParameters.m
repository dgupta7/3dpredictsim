% function res = f_WindLassSlackParameters(vars)
% vars = [0.14,0.2,1/500];
% k_PF = 8e+05;
% l_s = vars(1);
% q_s = vars(2);
% k_li = vars(3)*k_PF;
% Qs_mtp = [-30:15:30]*pi/180;
% 
% cWL = 0.025;
% a = 0.08207;
% b = 0.089638;
% phi0 = 2.493499;
% H0 = 0.027280;
% L0 = sqrt(a^2 + b^2 - 2*a*b*cos(phi0));
% 
% for i=1:length(Qs_mtp)
%     q_mtp = Qs_mtp(i);
%     l_0(i) = (1-cWL*(q_mtp*180/pi)/20)*L0; % foot arch length
%     h_0(i) = a*b/l_0(i)*sin(phi0);
%     q_mt_0(i) = acos( (a^2 + b^2 - l_0(i)^2)/(2*a*b) ) - phi0;
% 
%     dl_PF(i) = l_0(i) - l_s;
%     F_PFi = 11e3*(1-exp(-dl_PF(i)*80));
%     F_PF(i) = F_PFi*( tanh(F_PFi)+1 )/2;
%     
%     F_li(i) = k_li*(q_mt_0(i)-q_s)/h_0(i);
%     
% 
%     res(i) = F_PF(i) + F_li(i) - 5;
% end

cWL = 0.025;
R = vars(1);
L0 = 0.1628;
ls = vars(2);
l_toe = vars(3);

l_0_15 = (1-cWL*(15)/20)*L0;
l_0_30 = (1-cWL*(30)/20)*L0;
l_0_45 = (1-cWL*(45)/20)*L0;

l_15 = l_0_15 + R*15/180*pi + l_toe;
l_30 = l_0_30 + R*30/180*pi + l_toe;
l_45 = l_0_45 + R*45/180*pi + l_toe;


F_15 = 11e3*(1-exp((1-l_15/ls)*0.14/0.0125));
F_30 = 11e3*(1-exp((1-l_30/ls)*0.14/0.0125));
F_45 = 11e3*(1-exp((1-l_45/ls)*0.14/0.0125));
F_he = [F_15,F_30,F_45]/F_15;

F_15_1 = 7e5*(l_15-ls);
F_30_1 = 7e5*(l_30-ls);
F_45_1 = 7e5*(l_45-ls);
F_lin = [F_15_1,F_30_1,F_45_1]/F_15_1;

F_15_lit = 9.36*290/(l_15/ls);
F_30_lit = 11.7*290/(l_30/ls);
F_45_lit = 12.87*290/(l_45/ls);
F_lit = [F_15_lit,F_30_lit,F_45_lit]/F_15_lit;

dl_0 = ls/100*5;
F_15_tnh = 7e5*((l_15-ls) - dl_0*tanh((l_15-ls)/dl_0));
F_30_tnh = 7e5*((l_30-ls) - dl_0*tanh((l_30-ls)/dl_0));
F_45_tnh = 7e5*((l_45-ls) - dl_0*tanh((l_45-ls)/dl_0));
F_tnh = [F_15_tnh,F_30_tnh,F_45_tnh]/F_15_tnh;

x = [15,30,45];

figure
plot(x,F_he)
hold on
plot(x,F_lin)
plot(x,F_tnh)
plot(x,F_lit)
legend('hyperelastic','linear elastic',...
    'toe-in to lin el','literature')


res(1) = F_15/F_15_lit;
res(2) = F_30/F_30_lit;
res(3) = F_45/F_45_lit;


% end