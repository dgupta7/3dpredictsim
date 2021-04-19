

clear
clc


PF_stiffness = {'Natali2010'};
ls = 0.148;
f_PF_stiffness = f_getPlantarFasciaStiffnessModelCasADiFunction(PF_stiffness,'ls',ls);


q_mt = 0;
q_mtp = 0;

[M_mtj,M_mtp,M_PF,M_li,~,F,~,~,~,~,~] = getPassiveMtjMomentWindlass_v3(q_mt,0,q_mtp,f_PF_stiffness);

T_mtp = M_mtp - 10*q_mtp + 0.9254;

M_li_q0 = -M_PF;

%%
q_mt = 5*pi/180;
q_mtp = 10*pi/180;

[M_mtj,M_mtp,M_PF,M_li,~,F,~,~,~,~,~] = getPassiveMtjMomentWindlass_v3(q_mt,0,q_mtp,f_PF_stiffness);

T_mtj = -50;

M_li_q5 = T_mtj-M_PF;

%%
q_mt = -10*pi/180;
q_mtp = 30*pi/180;

[M_mtj,M_mtp,M_PF,M_li,~,F,~,~,~,~,~] = getPassiveMtjMomentWindlass_v3(q_mt,0,q_mtp,f_PF_stiffness);

T_mtj = 0;

M_li_q10 = -M_PF;

%%

q_mt = linspace(-15,7,1000)'*pi/180;

figure
plot([-10,0,5],[M_li_q10,M_li_q0,M_li_q5],'*')
hold on
grid on

%%

t1 = 5*pi/180;
c1 = 10;
c2 = 25;

t2 = 10*pi/180;
c3 = 2;
c4 = 40;
c5 = 1;

y = -c1*exp(c2*(q_mt-t1)) + c3*exp(-c4*(q_mt+t2)) + c5;

plot(q_mt*180/pi,y)












