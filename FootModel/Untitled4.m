clc
vars = [0.14,0,1/200];
[x,fval,~] = fsolve('f_WindLassSlackParameters',vars);

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
figure
plot(qt,la_r)


%%
M_li_0 = q_tmt_0(:,1)*kTMT_li;
F_li_0 = M_li_0./h0_fa(:,1);

figure
subplot(211)
plot(Qs_mtp*180/pi,M_li_0)
subplot(212)
plot(Qs_mtp*180/pi,F_li_0)

figure
plot(l0_fa(:,1),F_li_0)
