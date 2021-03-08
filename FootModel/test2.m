
clearvars -except 'S'
AddCasadiPaths();

%% Default settings
S = GetDefaultSettings(S);
clc
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
F  = external('F','Foot_v8.dll');
cd(pathmain);

%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
PathDefaultFunc = fullfile(pathCasADiFunctions,S.CasadiFunc_Folders);
f_PassiveMoments = Function.load(fullfile(PathDefaultFunc,'f_PassiveMoments'));
f_passiveWLTorques = Function.load(fullfile(PathDefaultFunc,'f_passiveWLTorques'));
f_forceEquilibrium_FtildeState_all_tendon = Function.load(fullfile(PathDefaultFunc,...
    'f_forceEquilibrium_FtildeState_all_tendon'));
f_lMT_vMT_dM = Function.load(fullfile(PathDefaultFunc,'f_lMT_vMT_dM'));
f_T12 = Function.load(fullfile(PathDefaultFunc,'f_T12'));

% Create a new function (since it doesn't exist for the full model)
% Function for 10 elements
ma_temp10 = SX.sym('ma_temp10',10);
ft_temp10 = SX.sym('ft_temp10',10);
J_sptemp10 = 0;
for i=1:length(ma_temp10)
    J_sptemp10 = J_sptemp10 + ma_temp10(i,1)*ft_temp10(i,1);
end
f_T10 = Function('f_T10',{ma_temp10,ft_temp10},{J_sptemp10});

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

%% Muscle information
% Muscles from one leg and from the back
muscleNames = {'glut_med1_r','glut_med2_r','glut_med3_r',...
    'glut_min1_r','glut_min2_r','glut_min3_r','semimem_r',...
    'semiten_r','bifemlh_r','bifemsh_r','sar_r','add_long_r',...
    'add_brev_r','add_mag1_r','add_mag2_r','add_mag3_r','tfl_r',...
    'pect_r','grac_r','glut_max1_r','glut_max2_r','glut_max3_r',......
    'iliacus_r','psoas_r','quad_fem_r','gem_r','peri_r',...
    'rect_fem_r','vas_med_r','vas_int_r','vas_lat_r','med_gas_r',...
    'lat_gas_r','soleus_r','tib_post_r','flex_dig_r','flex_hal_r',...
    'tib_ant_r','per_brev_r','per_long_r','per_tert_r','ext_dig_r',...
    'ext_hal_r','ercspn_r','intobl_r','extobl_r','ercspn_l',...
    'intobl_l','extobl_l'};
% Muscle indices for later use
pathmusclemodel = [pathRepo,'/MuscleModel'];
addpath(genpath(pathmusclemodel));
% Total number of muscles
NMuscle = length(muscleNames(1:end-3))*2;
pathpolynomial = fullfile(pathRepo,'Polynomials',S.PolyFolder); % default location 
tl = load([pathpolynomial,'/muscle_spanning_joint_INFO.mat']);
[~,mai] = MomentArmIndices(muscleNames(1:end-3),tl.muscle_spanning_joint_INFO);
% Muscles in foot
musif = find( ( tl.muscle_spanning_joint_INFO(:,5)==1 | tl.muscle_spanning_joint_INFO(:,6)==1 | ...
    tl.muscle_spanning_joint_INFO(:,7)==1 ) & tl.muscle_spanning_joint_INFO(:,4)==0 );
% musif = find(tl.muscle_spanning_joint_INFO(:,5)==1 | tl.muscle_spanning_joint_INFO(:,6)==1 | ...
%     tl.muscle_spanning_joint_INFO(:,7)==1);
NMf = length(musif);
for i=1:NMf
    muscleNamesFoot{i} = muscleNames{musif(i)};
end
tensions = getSpecificTensions(muscleNamesFoot);

%% Get Boundaries
load('D:\school\WTK\thesis\model\3dpredictsim\Results\MuscleModel\Pog_s1_bCst.mat','setup');
jointi = getJointi_tmt();

bounds_qs =  [[-30,0]*pi/180; % tibia rz
                    [0.2,0.6]; % tibia ty
                    [setup.bounds.Qs.lower(jointi.ankle.r), setup.bounds.Qs.upper(jointi.ankle.r)]...
                    *setup.scaling.Qs(jointi.ankle.r); % ankle
                    [setup.bounds.Qs.lower(jointi.subt.r), setup.bounds.Qs.upper(jointi.subt.r)]...
                    *setup.scaling.Qs(jointi.subt.r); % subt
                    [setup.bounds.Qs.lower(jointi.tmt.r), setup.bounds.Qs.upper(jointi.tmt.r)]...
                    *setup.scaling.Qs(jointi.tmt.r)]; % tmt

scale_qs = ones(size(bounds_qs,1),1);
for i=1:size(bounds_qs,1)
    scale_qs(i) = max(abs(bounds_qs(i,:)));
    bounds_qs(i,:) = bounds_qs(i,:)./scale_qs(i);
end

bounds_FTs = [zeros(NMf,1),ones(NMf,1)];
scale_FTs = 5*ones(NMf,1);

%% Build system to solve

% variables
Q_tib_rz = MX.sym('Q_tib_rz',1);
Q_tib_ty = MX.sym('Q_tib_ty',1);
Q_ankle = MX.sym('Q_ankle',1);
Q_subt = MX.sym('Q_subt',1);
Q_tmt = MX.sym('Q_tmt',1);
FT_tilde = MX.sym('FT_tilde',NMf,1);

% parameters
Q_mtp = MX.sym('Q_mtp',1);
F_tib_y = MX.sym('F_tib_y ',1);

% Get passive torques
Tau_pass_ankle = f_PassiveMoments(k_pass.ankle,theta.pass.ankle,Q_ankle,0);
Tau_pass_subt= f_PassiveMoments(k_pass.subt,theta.pass.subt,Q_subt,0);
Tau_pass_tmt = f_passiveWLTorques(Q_tmt,0,Q_mtp);

% Get muscle-tendon information
qin_r = MX.zeros(10,1); % adapt vector size to match full model
qin_r(1) = pi/2; % hip flex
qin_r(4) = -pi/2; % knee flex
qin_r(5) = Q_ankle;
qin_r(6) = Q_subt;
qin_r(7) = Q_mtp;
qdotin_r   = MX.zeros(10,1);
[lMTj_r,vMTj_r,MAj_r] = f_lMT_vMT_dM(qin_r,qdotin_r);
    
lMT_r = lMTj_r(musif);
vMT_r = vMTj_r(musif);

% Get muscle-tendon forces and derive Hill-equilibrium
akj = MX.zeros(92,1);  % adapt vector size to match full model
FTtildekj_nsc = MX.zeros(92,1);
dFTtildej_nsc = MX.zeros(92,1);
lMTj_lr = MX.zeros(92,1);
vMTj_lr = MX.zeros(92,1);
tensionsj_lr = MX.zeros(92,1);
for i=1:NMf
    FTtildekj_nsc((46+musif(i)),1) = FT_tilde(i)*scale_FTs(i);
    lMTj_lr((46+musif(i)),1) = lMT_r(i);
    vMTj_lr((46+musif(i)),1) = vMT_r(i);
    tensionsj_lr((46+musif(i)),1) = tensions(i);
end
[Hilldiffj,FTj,Fcej,Fpassj,Fisoj] =f_forceEquilibrium_FtildeState_all_tendon(akj(:,1),...
        FTtildekj_nsc(:,1),dFTtildej_nsc(:,1),lMTj_lr,vMTj_lr,tensionsj_lr);

Hilldiff = MX.zeros(NMf,1);
FT = MX.zeros(NMf,1);
for i=1:NMf
    Hilldiff(i) = Hilldiffj(46+musif(i));
    FT(i) = FTj(46+musif(i));
end

% Do not use mai 1 and 2, since they correspond with muscles also spanning
% the knee.
MAj.ankle.r      =  MAj_r(mai(5).mus.l(3:end)',5);
MAj.subt.r       =  MAj_r(mai(6).mus.l(3:end)',6);

% evaluate dynamics
qs = MX.zeros(nq,1);
qs(jointfi.tibia.rz) = Q_tib_rz*scale_qs(1);
qs(jointfi.tibia.ty) = Q_tib_ty*scale_qs(2);
qs(jointfi.ankle.r) = Q_ankle*scale_qs(3);
qs(jointfi.subt.r) = Q_subt*scale_qs(4);
qs(jointfi.tmt.r) = Q_tmt*scale_qs(5);
qs(jointfi.mtp.r) = Q_mtp;
qsqdots = MX.zeros(nq*2,1);
qsqdots(1:2:end,1) = qs;
A = MX.zeros(nq,1);

[Tj] = F([qsqdots(:,1);A(:,1)]);

% constraints
% Knee position
f1 = (Tj(jointfi.metatarsi_or(1),1) + Tj(jointfi.calcn_or(1),1))/2 - Tj(jointfi.tibia_or(1),1);
% Vertical force on knee
f2 = Tj(jointfi.tibia.ty,1) + F_tib_y;
% Ankle torque
Ft_ankle_r      = FT(:);
T_ankle_r       = f_T10(MAj.ankle.r,Ft_ankle_r);
f3 = Tj(jointfi.ankle.r,1) - (T_ankle_r + Tau_pass_ankle);
% Subtalar torque
Ft_subt_r       = FT(:);
T_subt_r        = f_T10(MAj.subt.r,Ft_subt_r);
f4 = Tj(jointfi.subt.r,1) - (T_subt_r + Tau_pass_subt);
% Tmt torque
f5 = Tj(jointfi.tmt.r,1) - Tau_pass_tmt;
% Hill difference
fh = Hilldiff;

ff = [f1;f2;f3;f4;f5;fh;];

f_foot = Function('f_foot',{[Q_tib_rz;Q_tib_ty;Q_ankle;Q_subt;Q_tmt;FT_tilde],Q_mtp,F_tib_y},{ff});

%%
% mtp angles to be considered
% Qs_mtp = [-30:5:45]*pi/180;
% Fs_tib = [0:100:1000];

% Qs_mtp = [-20,0,20]*pi/180;
Qs_mtp = [0]*pi/180;
Fs_tib = [0:100:500];

n_mtp = length(Qs_mtp);
n_tib = length(Fs_tib);


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
% positions
qs_opti = opti.variable(5,1);
opti.subject_to(bounds_qs(:,1) <= qs_opti(:,1));
opti.subject_to(qs_opti(:,1) <= bounds_qs(:,2));
% muscle-tendon forces
FTtilde_opti = opti.variable(NMf,1);
opti.subject_to(bounds_FTs(:,1) <= FTtilde_opti(:,1));
opti.subject_to(FTtilde_opti(:,1) <= bounds_FTs(:,2));
% parameters
qmtp = opti.parameter();
Ftib = opti.parameter();
% equality constraints
[constr] = f_foot([qs_opti;FTtilde_opti],qmtp,Ftib);
opti.subject_to(constr == 0);

% solver options
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.mu_strategy           = 'adaptive';
% options.ipopt.max_iter              = S.max_iter;
options.ipopt.linear_solver         = S.linear_solver;
options.ipopt.tol                   = 1*10^(-S.tol_ipopt);
opti.solver('ipopt', options);
     
temp = [];

%% run solver        
for i=1:n_mtp
    qs_init = [-10*pi/180;0.45;10*pi/180;0;0]./scale_qs;
    FTs_init = zeros(NMf,1);
    
    for j=1:n_tib
        %% solve the static situation
        % initial guess
        opti.set_initial(qs_opti,qs_init);
        opti.set_initial(FTtilde_opti,FTs_init);

        % set parameter values
        opti.set_value(qmtp, Qs_mtp(i));
        opti.set_value(Ftib, Fs_tib(j))
        % solve
        disp(['mtp: ' num2str(Qs_mtp(i)*180/pi) '° ; F tib: ' num2str(Fs_tib(j)) 'N']);
        sol = opti.solve();
        qs_sol = sol.value(qs_opti).*scale_qs;
        FTs_sol = sol.value(FTtilde_opti).*scale_FTs;
        
        % retrieve residuals
        f_res = full(f_foot([qs_sol./scale_qs;FTs_sol./scale_FTs],Qs_mtp(i),Fs_tib(j)));

        temp(:,end+1) = f_res;
        
        %% postprocess results
        % get states
        Qs(i,j,jointfi.tibia.rz) = qs_sol(1);
        Qs(i,j,jointfi.tibia.ty) = qs_sol(2);
        Qs(i,j,jointfi.ankle.r) = qs_sol(3);
        Qs(i,j,jointfi.subt.r) = qs_sol(4);
        Qs(i,j,jointfi.tmt.r) = qs_sol(5);
        Qs(i,j,jointfi.mtp.r) = Qs_mtp(i);

        QsQdots(i,j,1:2:end) = Qs(i,j,:);
        QsQdots(i,j,2:2:end) = Qdots(i,j,:);
        
        % call external function
        [T_res] = full(F([vertcat(squeeze(QsQdots(i,j,:)));vertcat(squeeze(Qddots(i,j,:)))]));
        
        % ground reaction forces
        GRF_calcn(i,j,:) = T_res(jointfi.calcn_GRF);
        GRF_metatarsi(i,j,:) = T_res(jointfi.metatarsi_GRF);

        % body origin positions
        toes_or(i,j,:) = T_res(jointfi.toes_or);
        metatarsi_or(i,j,:) = T_res(jointfi.metatarsi_or);
        calcn_or(i,j,:) = T_res(jointfi.calcn_or);
        talus_or(i,j,:) = T_res(jointfi.talus_or);
        tibia_or(i,j,:) = T_res(jointfi.tibia_or);
        
        % calculate arch length
        l_fa_ext(i,j) = norm(squeeze(toes_or(i,j,1:2)-calcn_or(i,j,1:2)));
        
        % calculate arch height
        x0 = metatarsi_or(i,j,1);
        x1 = calcn_or(i,j,1);
        x2 = toes_or(i,j,1);
        y0 = metatarsi_or(i,j,2);
        y1 = calcn_or(i,j,2);
        y2 = toes_or(i,j,2);
        h_fa_ext(i,j) = abs( (x2-x1)*(y1-y0) - (x1-x0)*(y2-y1) )/(sqrt( (x2-x1)^2 + (y2-y1)^2 ) );
        
        % call windlass function
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

%%

PlotResults_FootSim(R)






