
clear
clc

S.PF_stiffness = 'Natali2010';
S.R_mtth = 9.5e-3;
S.sf_PF = 10;                % multiply PF force with constant factor
S.PF_slack_length = 0.15; % (m) slack length
S.MT_li_nonl = 0;
S.kMT_li = 300;
S.kMTP = 5;

%%
import casadi.*
qin1     = SX.sym('qin_pass1',1);
qin2     = SX.sym('qin_pass2',1);

f_PF_stiffness = f_getPlantarFasciaStiffnessModelCasADiFunction(S.PF_stiffness,'ls',S.PF_slack_length);

[passWLTorques_mtj,passWLTorques_mtpj] = getPassiveMtjMomentWindlass_v3(qin1,0,qin2,f_PF_stiffness,S);

f_passiveWLTorques_mtj = Function('f_passiveWLTorques_mtj',{qin1,qin2}, ...
    {passWLTorques_mtj},{'qin1','qin2'},{'passWLTorques'});

if strcmp(S.PF_stiffness,'Cheng2008')
    offset = 1.11;
elseif strcmp(S.PF_stiffness,'Gefen2001')
    offset = 0.14;
elseif strcmp(S.PF_stiffness,'Ker1987')
    offset = 0.036;
elseif strcmp(S.PF_stiffness,'Natali2010')
    offset = 0.75;
elseif strcmp(S.PF_stiffness,'Song2011')
    offset = 0.037;
elseif strcmp(S.PF_stiffness,'linear')
    offset = 0.97;
elseif strcmp(S.PF_stiffness,'tanh')
    offset = 0.098;
else
    offset = 0;
end

k_pass.mtp = [-0.9 14.87 0.18 -70.08]';
theta.pass.mtp = [0 65/180*pi]';

Tau_pass = k_pass.mtp(1,1)*exp(k_pass.mtp(2,1)*(qin2-theta.pass.mtp(2,1))) + ...
    k_pass.mtp(3,1)*exp(k_pass.mtp(4,1)*(qin2-theta.pass.mtp(1,1)));

passTorques_mtpj = passWLTorques_mtpj - S.kMTP*qin2 + offset + Tau_pass;


f_passiveWLTorques_mtpj = Function('f_passiveWLTorques_mtpj',{qin1,qin2}, ...
    {passTorques_mtpj},{'qin1','qin2'},{'passWLTorques'});
 

%%  


Q_mtj = MX.sym('Q_mtj',1);
Q_mtp_gnd = MX.sym('Q_mtp_gnd',1);
Q_mtp = MX.sym('Q_mtp',1);
T_mtp_ext = MX.sym('T_mtp_ext',1);

T_mtj = f_passiveWLTorques_mtj(Q_mtj,Q_mtp);
T_mtp = f_passiveWLTorques_mtpj(Q_mtj,Q_mtp);

eq1 = -Q_mtp + Q_mtp_gnd - 0.4751*Q_mtj;
eq2 = T_mtj;
eq3 = T_mtp_ext - T_mtp;

f_opt_mtp = Function('f_opt',{[Q_mtj;Q_mtp;T_mtp_ext],Q_mtp_gnd},{[eq1;eq2;eq3]});



%% make solver
opti = casadi.Opti();
% positions
qs_opti = opti.variable(3,1);
opti.subject_to(-0.26<qs_opti(1)<0.26);
opti.subject_to(-0.8<qs_opti(2)<0.8);

% parameters
qmtp = opti.parameter();
% equality constraints
[constr] = f_opt_mtp(qs_opti,qmtp);
opti.subject_to(constr == 0);

opti.solver('ipopt');

opti.set_initial(qs_opti,[0,0,0]');


q_mtp = linspace(0,40,30);

for i=1:length(q_mtp)
    opti.set_value(qmtp, q_mtp(i)*pi/180);

    sol = opti.solve();
    qs_sol = sol.value(qs_opti);
    T(i) = full(qs_sol(3));
    q_mt(i) = full(qs_sol(1));
    q_mtpj(i) = full(qs_sol(2));
    
end

%%

figure
plot(q_mtp,T)

% figure
% plot(q_mtp,q_mt*180/pi)
% 
% figure
% plot(q_mtp,q_mtpj*180/pi)







