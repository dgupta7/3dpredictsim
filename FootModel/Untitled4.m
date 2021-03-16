clc

vars = [0.007, 0.16, 0.02];

% f_WindLassSlackParameters(vars);

% [x,fval,~] = fsolve('f_WindLassSlackParameters',vars);



%%
qt = [-45:5:45];
a = (qt+45)*pi/180;
R = 0.005;

lj = a*R;
lt = 0.005;
ls = 0.15;
la = ls - lj - lt;

la_0 = la(qt==0);
la_r = (la-la_0)/la_0;

% figure
% plot(qt,la_r)


%%
a = 0.08207;
b = 0.089638;
phi0 = 2.493499;
H0 = 0.027280;

L0 = sqrt(a^2 + b^2 - 2*a*b*cos(phi0));
cWL = 0.02;

for i=1:length(qt)
    l0_fa(i) = (1-cWL*(qt(i))/20)*L0; % foot arch length
    h0_fa(i) = a*b./l0_fa(i)*sin(phi0);
    q_mt_0(i) = acos( (a^2 + b^2 - l0_fa(i).^2)/(2*a*b) ) - phi0;

end
%%
kTMT_li = 150;
q0 = 20*pi/180;
M_li_0 = (q_mt_0(:)-q0)*kTMT_li;
F_li_0 = M_li_0./h0_fa(:);

figure
plot(qt,F_li_0)

figure
plot(l0_fa(:),F_li_0)

F_15 = interp1(qt,F_li_0,15);
F_30 = interp1(qt,F_li_0,30);
F_45 = interp1(qt,F_li_0,45);

F_15_lit = 9.36*290;
F_30_lit = 11.7*290;
F_45_lit = 12.87*290;

figure
plot([15,30,45],[F_15,F_30,F_45])
hold on
plot([15,30,45],[F_15_lit,F_30_lit,F_45_lit])



