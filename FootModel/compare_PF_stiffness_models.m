clear
clc

cWL = 0.020;
R = 0.02;
ls = 0.17;
l_toe = 0.015;

% linear stiffness
E = 350;
A0 = 290; % initial section (mm^2)
k1 = E*A0/ls;

% tanh stiffness parameter
k_1 = k1;
dl_0 = ls/100*1;

% quadratic stiffness
k2 = 5e7;

% polynomial stiffness parameters
a1 = -488737.9;
a2 = 2648898.5;
a3 = -5736967.6;
a4 = 6206986.7;
a5 = -3354935.1;
a6 = 724755.5;
Pr = 0.4; % Poisson ratio


% resulting stiffness parameters
kMT_li = 200;
q0_0 = 15*pi/180;
qf = 10;

% a = 0.08207;
% b = 0.089638;
% phi0 = 2.493499;
% H0 = 0.027280;

a = 0.0857;
b = 0.0932;
phi0 = 2.2799;
H0 = 0.0373;

L0 = sqrt(a^2 + b^2 - 2*a*b*cos(phi0));

qt = [-45:5:45];
ni = length(qt);

lz = linspace(1,1.13,1000);
l = lz*ls;
ni = length(l);

for i=1:ni
%     l0_fa(i) = (1-cWL*(qt(i))/20)*L0; % foot arch length
%     h0_fa(i) = a*b./l0_fa(i)*sin(phi0);
%     q_mt_0(i) = acos( (a^2 + b^2 - l0_fa(i).^2)/(2*a*b) ) - phi0;
%     l(i) = l0_fa(i) + R*qt(i)/180*pi + l_toe;
%     lz(i) = l(i)/ls;
    
    F(1,i) = 11e3*(1-exp((1-l(i)/ls)*0.14/0.0125));
    F_he(i) = F(1,i).*( tanh(F(1,i))+1 )/2;
    
    F(2,i) = k1*(l(i)-ls);
    F_lin(i) = F(2,i).*( tanh(F(2,i))+1 )/2;
    
    F(3,i) = k_1*((l(i)-ls) - dl_0*tanh((l(i)-ls)/dl_0));
    F_tnh(i) = F(3,i).*( tanh(F(3,i))+1 )/2;
    
    F(4,i) = k2*(l(i)-ls)^2;
    F_qdr(i) = F(4,i).*( tanh(F(4,i))+1 )/2;
    
    sgm(i) = a1*lz(i)^5 + a2*lz(i)^4 + a3*lz(i)^3 + a4*lz(i)^2 + a5*lz(i) + a6;
    A(i) = A0*(1-Pr*(lz(i)-1)).^2;
    F(5,i) = sgm(i)*A(i);
    F_ply(i) = F(5,i).*( tanh(F(5,i))+1 )/2;
    
    F(6,i) = 1e4*exp(l(i)-ls);
    F_exp(i) = F(6,i).*( tanh(F(6,i))+1 )/2;

%     q0(i) = q0_0 - qt(i)/qf *pi/180;
%     M_li(i) = (q_mt_0(i)-q0(i))*kMT_li;
%     F_li(i) = M_li(i)./h0_fa(i);

    
end

% figure
% plot(qt,F_lin/F_lin(10),'DisplayName','linear')
% hold on
% plot(qt,F_tnh/F_tnh(10),'DisplayName','stiffening')
% plot(qt,F_he/F_he(10),'DisplayName','hyperelastic')
% plot(qt,F_qdr/F_qdr(10),'DisplayName','quadratic')
% plot(qt,F_ply/F_ply(10),'DisplayName','5e order polynom')
% plot(qt,-F_li/abs(F_li(10)),'DisplayName','- other forces')
% legend('Location','best')
% xlabel('mtp angle(°)')
% ylabel('PF force (normalized at mtp = 0°)')


idx = find(lz>=1.02 & lz<=1.06);
k_lin = nanmean(F_lin(idx)./(l(idx)-ls));
k_tnh = nanmean(F_tnh(idx)./(l(idx)-ls));
k_he = nanmean(F_he(idx)./(l(idx)-ls));
k_qdr = nanmean(F_qdr(idx)./(l(idx)-ls));
k_ply = nanmean(F_ply(idx)./(l(idx)-ls));


figure
plot(lz,F_lin,'DisplayName','linear elastic')
hold on
plot(lz,F_tnh,'DisplayName','hypoelastic (tanh)')
plot(lz,F_he,'DisplayName','hyper elastic')
plot(lz,F_qdr,'DisplayName','hypoelastic (k \Deltal^2)')
plot(lz,F_ply,'DisplayName','hypoelastic (5e O)')
plot(lz,F_exp,'DisplayName','hypoelastic (exp)')
legend('Location','best')
xlabel('l/l_s (-)')
ylabel('PF force (N)')
% ylim([0,5000])

% figure
% hold on
% plot(qt,l/ls)





