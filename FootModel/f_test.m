function [res] = f_test(vars,f_PassiveMoments,f_passiveWLTorques,F,Q_mtp,F_tib_y,bounds_qs)


Q_tib_rz= vars(1);
Q_tib_ty= vars(2);
Q_ankle= vars(3);
Q_subt= vars(4);
Q_tmt= vars(5);


%% passive moments parameters
% copied from CreateCasADidiFunctions_all_tmt.m
k_pass.ankle = [-2.03 38.11 0.18 -12.12]';
theta.pass.ankle = [-0.74 0.52]';
k_pass.subt = [-60.21 16.32 60.21 -16.32]';
theta.pass.subt = [-0.65 0.65]';

%% Indices external function
% External function: F
% Joint torques.
jointfi.tibia.rz = 1;
jointfi.tibia.rx = 2;
jointfi.tibia.ry = 3;
jointfi.tibia.tx = 4;
jointfi.tibia.ty = 5;
jointfi.tibia.tz = 6;
jointfi.ankle.r = 7;
jointfi.subt.r = 8;
jointfi.tmt.r = 9;
jointfi.mtp.r = 10;
nq      = 10;
% Origin positions in ground frame
jointfi.tibia_or = 11:13;
jointfi.talus_or = 14:16;
jointfi.calcn_or = 17:19;
jointfi.metatarsi_or = 20:22;
jointfi.toes_or = 23:25;
% Ground reaction forces
jointfi.calcn_GRF = 26:28;
jointfi.metatarsi_GRF = 29:31;

% Get passive torques
Tau_pass_ankle = f_PassiveMoments(k_pass.ankle,theta.pass.ankle,Q_ankle,0);
Tau_pass_subt= f_PassiveMoments(k_pass.subt,theta.pass.subt,Q_subt,0);
Tau_pass_tmt = f_passiveWLTorques(Q_tmt,0,Q_mtp);


% evaluate dynamics
qs = zeros(nq,1);
qs(jointfi.tibia.rz) = Q_tib_rz;
qs(jointfi.tibia.ty) = Q_tib_ty;
qs(jointfi.ankle.r) = Q_ankle;
qs(jointfi.subt.r) = Q_subt;
qs(jointfi.tmt.r) = Q_tmt;
qs(jointfi.mtp.r) = Q_mtp;
qsqdots = zeros(nq*2,1);
qsqdots(1:2:end,1) = qs;
A = zeros(nq,1);
F0 = 0;

[Tj] = F([qsqdots(:,1);A(:,1);F0]);
% Tj = zeros(31,1);

% make function
% f1 = Tj(jointfi.metatarsi_or(1),1) - Tj(jointfi.tibia_or(1),1);
f2 = Tj(jointfi.tibia.ty,1) + F_tib_y; % vertical force on knee
f3 = Tj(jointfi.ankle.r,1) - Tau_pass_ankle;
f4 = Tj(jointfi.subt.r,1) - Tau_pass_subt;
f5 = Tj(jointfi.tmt.r,1) - Tau_pass_tmt;


f1 = (Tj(jointfi.metatarsi_or(1),1)+Tj(jointfi.calcn_or(1),1))/2 - Tj(jointfi.tibia_or(1),1); % position of knee
% f2 = Tj(jointfi.tibia.ty,1)+ F_tib_y; % vertical force on knee
% f3_1 = Tj(jointfi.ankle.r,1) - Tau_pass_ankle;
% f3_2 = Tj(jointfi.subt.r,1) - Tau_pass_subt;
% f3_3 = Tj(jointfi.tmt.r,1) - Tau_pass_tmt;
% f3 = f3_1^2 + f3_2^2 + f3_3^2;
% f4_1 = [Q_tib_rz;Q_tib_ty;Q_ankle;Q_subt;Q_tmt] - bounds_qs(:,1);
% f4_2 = (tanh(-f4_1*1e6)+1).*f4_1.^2;
% f4 = sum(f4_2);
% f5_1 = [Q_tib_rz;Q_tib_ty;Q_ankle;Q_subt;Q_tmt] - bounds_qs(:,2);
% f5_2 = (tanh(f5_1*1e6)+1).*f5_1.^2;
% f5 = sum(f5_2);

res = full([f1;f2;f3;f4;f5]);
end
