function [] = f_SimFoot2(S)

AddCasadiPaths();

%% Default settings
S = GetDefaultSettings(S);

%% Load external functions
import casadi.*
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB. More details
% about the external function can be found in the documentation.
pathmain        = pwd;
[pathRepo,~,~]  = fileparts(pathmain);
addpath(genpath(pathRepo));
% Loading external functions.
pathExternalFunctions = [pathRepo,'/ExternalFunctions'];
cd(pathExternalFunctions)
F  = external('F','Foot_v6.dll');
cd(pathmain);

%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
PathDefaultFunc = fullfile(pathCasADiFunctions,S.CasadiFunc_Folders);
f_PassiveMoments = Function.load(fullfile(PathDefaultFunc,'f_PassiveMoments'));
f_passiveWLTorques = Function.load(fullfile(PathDefaultFunc,'f_passiveWLTorques'));

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

%% Get Boundaries
load('D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst.mat','setup');
jointi = getJointi_tmt();

bounds_qs =  [[-30,0]*pi/180; % tibia rz
                    [0.2,0.6]; % tibia ty
                    [setup.bounds.Qs.lower(jointi.ankle.r), setup.bounds.Qs.upper(jointi.ankle.r)]*setup.scaling.Qs(jointi.ankle.r); % ankle
                    [setup.bounds.Qs.lower(jointi.subt.r), setup.bounds.Qs.upper(jointi.subt.r)]*setup.scaling.Qs(jointi.subt.r); % subt
                    [-15,15]*pi/180]; % tmt


%% Build system to solve

% variables
Q_tib_rz = MX.sym('Q_tib_rz',1);
Q_tib_ty = MX.sym('Q_tib_ty',1);
Q_ankle = MX.sym('Q_ankle',1);
Q_subt = MX.sym('Q_subt',1);
Q_tmt = MX.sym('Q_tmt',1);

% parameters
Q_mtp = MX.sym('Q_mtp',1);
F_tib_y = MX.sym('F_tib_y ',1);

% Get passive torques
Tau_pass_ankle = f_PassiveMoments(k_pass.ankle,theta.pass.ankle,Q_ankle,0);
Tau_pass_subt= f_PassiveMoments(k_pass.subt,theta.pass.subt,Q_subt,0);
Tau_pass_tmt = f_passiveWLTorques(Q_tmt,0,Q_mtp);

% evaluate dynamics
qs = MX.sym('qs',nq);
qs(:) = 0;
qs(jointfi.tibia.rz) = Q_tib_rz;
qs(jointfi.tibia.ty) = Q_tib_ty;
qs(jointfi.ankle.r) = Q_ankle;
qs(jointfi.subt.r) = Q_subt;
qs(jointfi.tmt.r) = Q_tmt;
qs(jointfi.mtp.r) = Q_mtp;
qsqdots = MX.sym('qsqdots',nq*2);
qsqdots(:) = 0;
qsqdots(1:2:end,:) = qs;
A = MX.sym('qdds',nq);
A(:) = 0;
F0 = MX.sym('F0',1);
F0(:) = 0;

[Tj] = F([qsqdots(:,1);A(:,1);F0]);
% Tj = zeros(31,1);

% positions
% toes_or_sx = Tj(jointi.toes_or);
% metatarsi_or_sx = Tj(jointfi.metatarsi_or);
% calcn_or_sx = Tj(jointfi.calcn_or);
% talus_or_sx = Tj(jointi.talus_or);
% tibia_or_sx = Tj(jointfi.tibia_or);

% make function
% f1 = metatarsi_or_sx(1) - tibia_or_sx(1); % position of knee
f1 = Tj(jointfi.metatarsi_or) - Tj(jointfi.tibia_or);
f2 = Tj(jointfi.tibia.ty,1) + F_tib_y; % vertical force on knee
f3 = Tj(jointfi.ankle.r,1) - Tau_pass_ankle;
f4 = Tj(jointfi.subt.r,1) - Tau_pass_subt;
f5 = Tj(jointfi.tmt.r,1) - Tau_pass_tmt;

% f1 = 0;
% f2 = 0;
% f3 = 0;
% f4 = 0;
% f5 = 0;

% ff2 = f1^2 + f2^2 + f3^2 + f4^2 + f5^2;

f_foot = Function('f_foot',{[Q_tib_rz;Q_tib_ty;Q_ankle;Q_subt;Q_tmt],Q_mtp,F_tib_y},{[f1;f2;f3;f4;f5]});
% f_foot = Function('f_foot',{[Q_tib_rz;Q_tib_ty;Q_ankle;Q_subt;Q_tmt],Q_mtp,F_tib_y},{ff2});

% f_foot_sol = rootfinder('f_foot_sol','newton',f_foot);
% disp(f_foot_sol)


%%
% mtp angles to be considered
% Qs_mtp = [-30:5:45]*pi/180;
% Fs_tib = [0:100:1000];

% Qs_mtp = [-20,0,20]*pi/180;
Qs_mtp = [0]*pi/180;
Fs_tib = [0];

n_mtp = length(Qs_mtp);
n_tib = length(Fs_tib);

% % initial guess for solver
% 
% optim_options = optimset('Display','off');

% Windlass parameters
kTMT_li = 1.5/(pi/180)/5;
kTMT_PF = S.kTMT;
dTMT = S.dTMT;
cWL = S.cWL;

% Declare arrays for postprocessing results
Qs = zeros(n_mtp,n_tib,nq);
Qdots = zeros(n_mtp,n_tib,nq);
Qddots = zeros(n_mtp,n_tib,nq);
QsQdots = zeros(n_mtp,n_tib,nq*2);
GRF_calcn = zeros(n_mtp,n_tib,3);
GRF_metatarsi = zeros(n_mtp,n_tib,3);
toes_or = zeros(n_mtp,n_tib,3);
calcn_or = zeros(n_mtp,n_tib,3);
metatarsi_or = zeros(n_mtp,n_tib,3);
talus_or = zeros(n_mtp,n_tib,3);
tibia_or = zeros(n_mtp,n_tib,3);
l_fa_ext = zeros(n_mtp,n_tib);
h_fa_ext = zeros(n_mtp,n_tib);
M = zeros(n_mtp,n_tib);
M_PF = zeros(n_mtp,n_tib);
F_PF = zeros(n_mtp,n_tib);
l_fa = zeros(n_mtp,n_tib);
h_fa = zeros(n_mtp,n_tib);
l0_fa = zeros(n_mtp,n_tib);
h0_fa = zeros(n_mtp,n_tib);
q_tmt_0 = zeros(n_mtp,n_tib);

%% make solver
opti = casadi.Opti();
qs_opti = opti.variable(5,1);
% opti.subject_to(bounds_qs(:,1) <= qs_opti(:,1));
% opti.subject_to(qs_opti(:,1) <= bounds_qs(:,2));
% qmtp = opti.parameter();
% Ftib = opti.parameter();
qmtp = 0;
Ftib = 0;
[fsi] = f_foot(qs_opti,qmtp,Ftib);
opti.subject_to(fsi == 0);
% opti.minimize(fsi(1,1)^2 + fsi(2,1)^2 + fsi(3,1)^2 + fsi(4,1)^2 + fsi(5,1)^2);
% opti.minimize(fsi);
opti.solver('ipopt');

        
        
for i=1:n_mtp
%     vars_init = zeros(5,1);
%     vars_init(2) = 0.451;
    vars_init = [-10*pi/180;0.45;10*pi/180;0;0];
    for j=1:n_tib
        %% solve the static situation
        opti.set_initial(qs_opti,vars_init);
        opti.set_value(qmtp, Qs_mtp(i));
        opti.set_value(Ftib, Fs_tib(j))
        
        sol = opti.solve();
        x = sol.value(qs_opti);

%     prob = struct('f', opti.f, 'x', opti.x, 'g',opti.g );
%     solver = nlpsol('solver', 'ipopt',prob);
%     lbx = bounds_qs(:,1);
%     ubx = bounds_qs(:,2);
%     lbg = zeros(5,1);
%     ubg = zeros(5,1);
%     sol = solver('x0',guess,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg); 
%     x = full(sol.x);
    
%         vars_init = x; % use solution as ig for next (higher) force
        
        %% postprocess results
        % get states
        Qs(i,j,jointfi.tibia.rz) = x(1);
        Qs(i,j,jointfi.tibia.ty) = x(2);
        Qs(i,j,jointfi.ankle.r) = x(3);
        Qs(i,j,jointfi.subt.r) = x(4);
        Qs(i,j,jointfi.tmt.r) = x(5);
        Qs(i,j,jointfi.mtp.r) = Qs_mtp(i);

        QsQdots(i,j,1:2:end) = Qs(i,j,:);
        QsQdots(i,j,2:2:end) = Qdots(i,j,:);
        
        % call external function
        [T_res] = full(F([vertcat(squeeze(QsQdots(i,j,:)));vertcat(squeeze(Qddots(i,j,:)));0]));
        
        % get GRF
        GRF_calcn(i,j,:) = T_res(jointfi.calcn_GRF);
        GRF_metatarsi(i,j,:) = T_res(jointfi.metatarsi_GRF);

        toes_or(i,j,:) = T_res(jointfi.toes_or);
        metatarsi_or(i,j,:) = T_res(jointfi.metatarsi_or);
        calcn_or(i,j,:) = T_res(jointfi.calcn_or);
        talus_or(i,j,:) = T_res(jointfi.talus_or);
        tibia_or(i,j,:) = T_res(jointfi.tibia_or);
        
        l_fa_ext(i,j) = norm(squeeze(toes_or(i,j,1:2)-calcn_or(i,j,1:2)));
        h_fa_ext(i,j) = metatarsi_or(i,j,2)-calcn_or(i,j,2);
        
        [Mi, M_PFi,F_PFi,~,~,li,l0i,L0,hi,h0i,H0,q_tmt_0i] = ...
                getPassiveTmtjMomentWindlass(Qs(i,j,jointfi.tmt.r),0,Qs(i,j,jointfi.mtp.r),...
                kTMT_li,kTMT_PF,dTMT,S.subject,cWL);

        M(i,j) = Mi;
        M_PF(i,j) = M_PFi;
        l_fa(i,j) = li;
        h_fa(i,j) = hi;
        F_PF(i,j) = F_PFi;
        l0_fa(i,j) = l0i;
        h0_fa(i,j) = h0i;
        q_tmt_0(i,j) = q_tmt_0i;

    
    end
end


%% save results
R.S = S;
R.Qs_mtp = Qs_mtp;
R.Qs = Qs;
R.jointfi = jointfi;
R.Fs_tib = Fs_tib;
R.GRF_calcn = GRF_calcn;
R.GRF_metatarsi = GRF_metatarsi;
R.M_WL = M;
R.M_PF = M_PF;
R.F_PF = F_PF;
R.l_fa = l_fa;
R.l_fa_ext = l_fa_ext;
R.l0_fa = l0_fa;
R.L0 = L0;
R.h_fa = h_fa;
R.h_fa_ext = h_fa_ext;
R.h0_fa = h0_fa;
R.H0 = H0;
R.Q_tmt_0 = q_tmt_0;
R.toes_or = toes_or;
R.metatarsi_or = metatarsi_or;
R.calcn_or = calcn_or;
R.talus_or = talus_or;
R.tibia_or = tibia_or;

% OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
% FilenameAnalysis = fullfile(OutFolder,[S.savename '_pp.mat']);
% save(FilenameAnalysis,'R');



PlotResults_FootSim(R)






end