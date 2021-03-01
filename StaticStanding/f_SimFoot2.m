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

%% Indices external function
% External function: F
% Joint torques.
jointfi.tibia.rx = 1;
jointfi.tibia.ry = 2;
jointfi.tibia.rz = 3;
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

%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
PathDefaultFunc = fullfile(pathCasADiFunctions,S.CasadiFunc_Folders);
f_AllPassiveTorques = Function.load(fullfile(PathDefaultFunc,'f_AllPassiveTorques'));

%% Build simplified passive torque function
% Indices of passive torques casadi function
jointi = getJointi_tmt();
residualsi = jointi.pelvis.tilt:jointi.elb.r; % all
nqall = length(residualsi); % all

Qs_p = SX.sym('QS_pass',4);
Qs_nsc = SX.sym('Qs_nsc',nqall);
Qs_nsc(:) = 0;
Qdots_nsc = zeros(nqall,1);

Qs_nsc(jointi.ankle.r) = Qs_p(1);
Qs_nsc(jointi.subt.r) = Qs_p(2);
Qs_nsc(jointi.tmt.r) = Qs_p(3);
Qs_nsc(jointi.mtp.r) = Qs_p(4);

% Get passive joint torques
Tau_passj_all = full(f_AllPassiveTorques(Qs_nsc(:,1),Qdots_nsc(:,1)));
Tau_p = [Tau_passj_all(10),Tau_passj_all(12),Tau_passj_all(14)];


f_getPassT_Foot = Function('f_getPassT_Foot',{Qs_p},{Tau_p});



%% Build system to solve

% variables
vars = MX.sym('vars',5);
Q_ankle = vars(1);
Q_subt = vars(2);
Q_tmt = vars(3);
Q_tib_rz = vars(4);
Q_tib_ty = vars(5);

% Q_ankle = MX.sym('Q',1);
% Q_subt = MX.sym('Q',1);
% Q_tmt = MX.sym('Q',1);
% Q_tib_rz = MX.sym('Q',1);
% Q_tib_ty = MX.sym('Q',1);

% inputs
Q_mtp = MX.sym('Q_mtp',1);
F_tib_y = MX.sym('F_tib_y ',1);

% Get passive torques from previously made function
Qs_pass = [Q_ankle,Q_subt,Q_tmt,Q_mtp]; 
Tau_pass = f_getPassT_Foot(Qs_pass);

% evaluate dynamics
qs = MX.sym('qs',nq);
qs(:) = 0;
qs(jointfi.tibia.rz) = Q_tib_rz;
qs(jointfi.tibia.ty) = Q_tib_ty;
qs(jointfi.ankle.r) = Q_ankle;
qs(jointfi.subt.r) = Q_subt;
qs(jointfi.tmt.r) = Q_tmt;
qs(jointfi.mtp.r) = Q_mtp;
A = zeros(nq,1);
qsqdots = MX.sym('qsqdots',nq*2);
qsqdots(:) = 0;
qsqdots(1:2:end,:) = qs;

[Tj] = F([qsqdots(:,1);A(:,1);0]);

% positions
% toes_or_sx = Tj(jointi.toes_or);
metatarsi_or_sx = Tj(jointfi.metatarsi_or);
calcn_or_sx = Tj(jointfi.calcn_or);
% talus_or_sx = Tj(jointi.talus_or);
tibia_or_sx = Tj(jointfi.tibia_or);

% make function
f1 = (metatarsi_or_sx(2)+calcn_or_sx(2))/2 - tibia_or_sx(2); % position of knee
f2 = Tj(jointfi.tibia.ty,1)+ F_tib_y; % vertical force on knee
f3 = Tj(jointfi.ankle.r,1) - Tau_pass(1); % ankle
f4 = Tj(jointfi.subt.r,1) - Tau_pass(2); % subt
f5 = Tj(jointfi.tmt.r,1) - Tau_pass(3); % tmt

fs = Function('fs',{vars,Q_mtp,F_tib_y},{[f1;f2;f3;f4;f5]});


% G = rootfinder('G','newton',fs);
% disp(G)


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

%%

for i=1:n_mtp
%     vars_init = zeros(nq-1,1);
%     vars_init(jointfi.tibia.ty) = 0.4;
    vars_init = [0;0;0;0;0.4];
    for j=1:n_tib
        %% solve the static situation
%         [x,fval,~]=fsolve('f_Foot',vars_init,optim_options,F,PassT_Foot,jointfi,Qs_mtp(i),Fs_tib(j));
%         [x] = G(vars_init,Qs_mtp(i),Fs_tib(j));

        x_MX = MX.sym('x',5,1);
        [fsi] = fs(x_MX,Qs_mtp(i),Fs_tib(j));
        cost = fsi(1)^2 + fsi(2)^2 + fsi(3)^2 + fsi(4)^2 + fsi(5)^2;
        constr = fsi(1);
        nlp = struct('x', x_MX, 'f', cost, 'g', constr);
        solver = nlpsol('solver', 'ipopt', nlp);
        res = solver('x0' , vars_init, 'lbg',0,'ubg',0);
        x = full(res);

        vars_init = x; %use solution as ig for next (higher) force
        %% postprocess results
        % get states
        Qs(i,j,1:9) = x;
        Qs(i,j,10) = Qs_mtp(i);

        QsQdots(i,j,1:2:end) = Qs(i,j,:);
        QsQdots(i,j,2:2:end) = Qdots(i,j,:);

        % call external function
        [T_res] = full(F([vertcat(squeeze(QsQdots(i,j,:)));vertcat(squeeze(Qddots(i,j,:)));Fs_tib(j)]));
        
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

OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
FilenameAnalysis = fullfile(OutFolder,[S.savename '_pp.mat']);
save(FilenameAnalysis,'R');


%%
PlotResults_FootSim(R)






end