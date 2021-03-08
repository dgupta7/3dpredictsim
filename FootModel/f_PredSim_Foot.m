function [] = f_PredSim_Foot(S)

%% Adding the casadi path seems to be needed to run processes in batch
AddCasadiPaths();

%% Default settings
S = GetDefaultSettings(S);

%% User inputs (typical settings structure)
% settings for optimization
N           = S.N;          % number of mesh intervals

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
setup.derivatives =  'AD'; % Algorithmic differentiation
pathExternalFunctions = [pathRepo,'/ExternalFunctions'];
cd(pathExternalFunctions)
F  = external('F','Foot_v8.dll');
cd(pathmain);

%% Indices external function
% Indices of the elements in the external functions
% Indices for foot model
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
nq.foot      = 10;
% Origin positions in ground frame
jointfi.tibia_or = 11:13;
jointfi.talus_or = 14:16;
jointfi.calcn_or = 17:19;
jointfi.metatarsi_or = 20:22;
jointfi.toes_or = 23:25;
% Ground reaction forces
jointfi.calcn_GRF = 26:28;
jointfi.metatarsi_GRF = 29:31;

% Indices of full model
jointi = getJointi_tmt();

% Vectors of indices for later use
residualsi          = jointi.pelvis.tilt:jointi.elb.r; % all
ground_pelvisi      = jointi.pelvis.tilt:jointi.pelvis.tz; % ground-pelvis
trunki              = jointi.trunk.ext:jointi.trunk.rot; % trunk
armsi               = jointi.sh_flex.l:jointi.elb.r; % arms
mtpi                = jointi.mtp.l:jointi.mtp.r; % mtps
residuals_noarmsi   = jointi.pelvis.tilt:jointi.trunk.rot; % all but arms
% Number of degrees of freedom for later use
nq.all      = length(residualsi); % all
nq.abs      = length(ground_pelvisi); % ground-pelvis
nq.trunk    = length(trunki); % trunk
nq.arms     = length(armsi); % arms
nq.mtp      = length(mtpi); % arms
nq.leg      = 10; % #joints needed for polynomials


%% Collocation scheme
% We use a pseudospectral direct collocation method, i.e. we use Lagrange
% polynomials to approximate the state derivatives at the collocation
% points in each mesh interval. We use d=3 collocation points per mesh
% interval and Radau collocation points.
pathCollocationScheme = [pathRepo,'/CollocationScheme'];
addpath(genpath(pathCollocationScheme));
d = 3; % degree of interpolating polynomial
method = 'radau'; % collocation method
[~,C,D,B] = CollocationScheme(d,method);

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

% Parameters for activation dynamics
tact = 0.015; % Activation time constant
tdeact = 0.06; % Deactivation time constant

%% Metabolic energy model parameters
% We extract the specific tensions and slow twitch rations.
pathMetabolicEnergy = [pathRepo,'/MetabolicEnergy'];
addpath(genpath(pathMetabolicEnergy));
tension = getSpecificTensions(muscleNames(1:end-3));
tensions = [tension;tension];
pctst = getSlowTwitchRatios(muscleNames(1:end-3));
pctsts = [pctst;pctst];

%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
PathDefaultFunc = fullfile(pathCasADiFunctions,S.CasadiFunc_Folders);
f_ArmActivationDynamics = Function.load(fullfile(PathDefaultFunc,'f_ArmActivationDynamics'));
f_FiberLength_TendonForce_tendon = Function.load(fullfile(PathDefaultFunc,'f_FiberLength_TendonForce_tendon'));
f_FiberVelocity_TendonForce_tendon = Function.load(fullfile(PathDefaultFunc,'f_FiberVelocity_TendonForce_tendon'));
f_forceEquilibrium_FtildeState_all_tendon = Function.load(fullfile(PathDefaultFunc,'f_forceEquilibrium_FtildeState_all_tendon'));
f_J2    = Function.load(fullfile(PathDefaultFunc,'f_J2'));
f_J25   = Function.load(fullfile(PathDefaultFunc,'f_J25'));
f_J23   = Function.load(fullfile(PathDefaultFunc,'f_J23'));
f_J8    = Function.load(fullfile(PathDefaultFunc,'f_J8'));
f_J92   = Function.load(fullfile(PathDefaultFunc,'f_J92'));
f_J92exp = Function.load(fullfile(PathDefaultFunc,'f_J92exp'));
f_Jnn2  = Function.load(fullfile(PathDefaultFunc,'f_Jnn2'));
f_lMT_vMT_dM = Function.load(fullfile(PathDefaultFunc,'f_lMT_vMT_dM'));
f_MtpActivationDynamics = Function.load(fullfile(PathDefaultFunc,'f_MtpActivationDynamics'));
f_T12 = Function.load(fullfile(PathDefaultFunc,'f_T12'));
f_T13 = Function.load(fullfile(PathDefaultFunc,'f_T13'));
f_T27 = Function.load(fullfile(PathDefaultFunc,'f_T27'));
f_T6 = Function.load(fullfile(PathDefaultFunc,'f_T6'));
f_AllPassiveTorques = Function.load(fullfile(PathDefaultFunc,'f_AllPassiveTorques'));
fgetMetabolicEnergySmooth2004all = Function.load(fullfile(PathDefaultFunc,'fgetMetabolicEnergySmooth2004all'));

% file with mass of muscles
MassFile = fullfile(PathDefaultFunc,'MassM.mat');
if exist(MassFile,'file')
    MuscleMass = load(MassFile);
else
    MassFile = fullfile(pathCasADiFunctions,'MassM.mat'); % default muscle mass
    MuscleMass =load(MassFile);
end



%% Get bounds and initial guess

% Kinematics file for bounds -- input arguments
IKfile_bounds = fullfile(pathRepo, S.IKfile_Bounds);

% We extract experimental data to set bounds and initial guesses if needed
joints = {'pelvis_tilt','pelvis_list','pelvis_rotation','pelvis_tx',...
    'pelvis_ty','pelvis_tz','hip_flexion_l','hip_adduction_l',...
    'hip_rotation_l','hip_flexion_r','hip_adduction_r','hip_rotation_r',...
    'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
    'subtalar_angle_l','subtalar_angle_r','tmt_angle_l','tmt_angle_r',...
    'mtp_angle_l','mtp_angle_r',...
    'lumbar_extension','lumbar_bending','lumbar_rotation','arm_flex_l',...
    'arm_add_l','arm_rot_l','arm_flex_r','arm_add_r','arm_rot_r',...
    'elbow_flex_l','elbow_flex_r'};

Qs_walk          = getIK(IKfile_bounds,joints);
[bounds,scaling] = getBounds_all_tmt(Qs_walk,NMuscle,nq,jointi,S.v_tgt);

% adapt bounds based on user input
bounds = AdaptBounds(bounds,S,mai);

bounds_qs_foot =   [[-30,0]*pi/180; % tibia rz
                    [0,0]; % tibia rx
                    [0,0]; % tibia ry
                    [0,0]; % tibia tx
                    [0.2,0.6]; % tibia ty
                    [0,0]; % tibia tz
                    [bounds.Qs.lower(jointi.ankle.r), bounds.Qs.upper(jointi.ankle.r)]*scaling.Qs(jointi.ankle.r); % ankle
                    [bounds.Qs.lower(jointi.subt.r), bounds.Qs.upper(jointi.subt.r)]*scaling.Qs(jointi.subt.r); % subt
                    [bounds.Qs.lower(jointi.tmt.r), bounds.Qs.upper(jointi.tmt.r)]*scaling.Qs(jointi.tmt.r); % tmt
                    [bounds.Qs.lower(jointi.mtp.r), bounds.Qs.upper(jointi.mtp.r)]*scaling.Qs(jointi.mtp.r)]; % mtp

bounds_qdots_foot =   [[-30,30]*pi/180; % tibia rz
                    [0,0]; % tibia rx
                    [0,0]; % tibia ry
                    [0,0]; % tibia tx
                    [-1,1]; % tibia ty
                    [0,0]; % tibia tz
                    [bounds.Qdots.lower(jointi.ankle.r), bounds.Qdots.upper(jointi.ankle.r)]*scaling.Qdots(jointi.ankle.r); % ankle
                    [bounds.Qdots.lower(jointi.subt.r), bounds.Qdots.upper(jointi.subt.r)]*scaling.Qdots(jointi.subt.r); % subt
                    [bounds.Qdots.lower(jointi.tmt.r), bounds.Qdots.upper(jointi.tmt.r)]*scaling.Qdots(jointi.tmt.r); % tmt
                    [bounds.Qdots.lower(jointi.mtp.r), bounds.Qdots.upper(jointi.mtp.r)]*scaling.Qdots(jointi.mtp.r)]; % mtp

bounds_qddots_foot =   [[-300,300]*pi/180; % tibia rz
                    [0,0]; % tibia rx
                    [0,0]; % tibia ry
                    [0,0]; % tibia tx
                    [-10,10]; % tibia ty
                    [0,0]; % tibia tz
                    [bounds.Qdotdots.lower(jointi.ankle.r), bounds.Qdotdots.upper(jointi.ankle.r)]*scaling.Qdotdots(jointi.ankle.r); % ankle
                    [bounds.Qdotdots.lower(jointi.subt.r), bounds.Qdotdots.upper(jointi.subt.r)]*scaling.Qdotdots(jointi.subt.r); % subt
                    [bounds.Qdotdots.lower(jointi.tmt.r), bounds.Qdotdots.upper(jointi.tmt.r)]*scaling.Qdotdots(jointi.tmt.r); % tmt
                    [bounds.Qdotdots.lower(jointi.mtp.r), bounds.Qdotdots.upper(jointi.mtp.r)]*scaling.Qdotdots(jointi.mtp.r)]; % mtp


% The initial guess depends on the settings
pathIG = [pathRepo,'/IG'];
addpath(genpath(pathIG));

guess = getGuess_QR_opti_int_tmt(N,nq,NMuscle,scaling,S.v_tgt,jointi,d,S.IG_PelvisY);

% adapt guess so that it fits within the bounds
guess = AdaptGuess_UserInput(guess,bounds,S);

%% external force on tibia
% This force is defined as going from 0 to full amplitude to 0 over one
% periode. A static offset can be added.
t_F_tib_y = linspace(0,S.FootModel.F_T,N+1);
F_tib_y = S.FootModel.F_Offset + S.FootModel.F_Amplitude/2*(1-cos(2*pi*t_F_tib_y/S.FootModel.F_T));


%% Index helpers
% get help indexes for left and right leg and for symmetry constraint
[IndexLeft,IndexRight,QsInvA,QsInvB,QdotsInvA,...
    QdotsInvB,orderQsOpp] = GetIndexHelper_tmt(S,jointi);

%% OCP create variables and bounds
% using opti
opti = casadi.Opti();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define static parameters
% Final time
tf = S.FootModel.F_T;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define states
% Muscle activations at mesh points
a = opti.variable(NMuscle,N+1);
opti.subject_to(bounds.a.lower'*ones(1,N+1) < a < ...
    bounds.a.upper'*ones(1,N+1));
opti.set_initial(a, guess.a');
% Muscle activations at collocation points
a_col = opti.variable(NMuscle,d*N);
opti.subject_to(bounds.a.lower'*ones(1,d*N) < a_col < ...
    bounds.a.upper'*ones(1,d*N));
opti.set_initial(a_col, guess.a_col');
% Muscle-tendon forces at mesh points
FTtilde = opti.variable(NMuscle,N+1);
opti.subject_to(bounds.FTtilde.lower'*ones(1,N+1) < FTtilde < ...
    bounds.FTtilde.upper'*ones(1,N+1));
opti.set_initial(FTtilde, guess.FTtilde');
% Muscle-tendon forces at collocation points
FTtilde_col = opti.variable(NMuscle,d*N);
opti.subject_to(bounds.FTtilde.lower'*ones(1,d*N) < FTtilde_col < ...
    bounds.FTtilde.upper'*ones(1,d*N));
opti.set_initial(FTtilde_col, guess.FTtilde_col');
% Qs at mesh points
Qs = opti.variable(nq.foot,N+1);

% Qs at collocation points
Qs_col = opti.variable(nq.foot,d*N);
opti.subject_to(bounds.QsQdots.lower(1:2:end)'*ones(1,d*N) < Qs_col < ...
    bounds.QsQdots.upper(1:2:end)'*ones(1,d*N));
opti.set_initial(Qs_col, guess.QsQdots_col(:,1:2:end)');
% Qdots at mesh points
Qdots = opti.variable(nq.all,N+1);
opti.subject_to(bounds.QsQdots.lower(2:2:end)'*ones(1,N+1) < Qdots < ...
    bounds.QsQdots.upper(2:2:end)'*ones(1,N+1));
opti.set_initial(Qdots, guess.QsQdots(:,2:2:end)');
% Qdots at collocation points
Qdots_col = opti.variable(nq.all,d*N);
opti.subject_to(bounds.QsQdots.lower(2:2:end)'*ones(1,d*N) < Qdots_col < ...
    bounds.QsQdots.upper(2:2:end)'*ones(1,d*N));
opti.set_initial(Qdots_col, guess.QsQdots_col(:,2:2:end)');
% Arm activations at mesh points
a_a = opti.variable(nq.arms,N+1);
opti.subject_to(bounds.a_a.lower'*ones(1,N+1) < a_a < ...
    bounds.a_a.upper'*ones(1,N+1));
opti.set_initial(a_a, guess.a_a');
% Arm activations at collocation points
a_a_col = opti.variable(nq.arms,d*N);
opti.subject_to(bounds.a_a.lower'*ones(1,d*N) < a_a_col < ...
    bounds.a_a.upper'*ones(1,d*N));
opti.set_initial(a_a_col, guess.a_a_col');
% Arm activations at mesh points
a_mtp = opti.variable(nq.mtp,N+1);
opti.subject_to(bounds.a_mtp.lower'*ones(1,N+1) < a_mtp < ...
    bounds.a_mtp.upper'*ones(1,N+1));
opti.set_initial(a_mtp, guess.a_mtp');
% Mtp activations at collocation points
a_mtp_col = opti.variable(nq.mtp,d*N);
opti.subject_to(bounds.a_mtp.lower'*ones(1,d*N) < a_mtp_col < ...
    bounds.a_mtp.upper'*ones(1,d*N));
opti.set_initial(a_mtp_col, guess.a_mtp_col');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define controls
% Time derivative of muscle activations (states) at mesh points
vA = opti.variable(NMuscle, N);
opti.subject_to(bounds.vA.lower'*ones(1,N) < vA < ...
    bounds.vA.upper'*ones(1,N));
opti.set_initial(vA, guess.vA');
% Arm excitations
e_a = opti.variable(nq.arms, N);
opti.subject_to(bounds.e_a.lower'*ones(1,N) < e_a < ...
    bounds.e_a.upper'*ones(1,N));
opti.set_initial(e_a, guess.e_a');
% Mtp excitations
e_mtp = opti.variable(nq.mtp, N);
opti.subject_to(bounds.e_mtp.lower'*ones(1,N) < e_mtp < ...
    bounds.e_mtp.upper'*ones(1,N));
opti.set_initial(e_mtp, guess.e_mtp');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define "slack" controls
% Time derivative of muscle-tendon forces (states) at collocation points
dFTtilde_col = opti.variable(NMuscle, d*N);
opti.subject_to(bounds.dFTtilde.lower'*ones(1,d*N) < dFTtilde_col < ...
    bounds.dFTtilde.upper'*ones(1,d*N));
opti.set_initial(dFTtilde_col, guess.dFTtilde_col');
% Time derivative of Qdots (states) at collocation points
A_col = opti.variable(nq.all, d*N);
opti.subject_to(bounds.Qdotdots.lower'*ones(1,d*N) < A_col < ...
    bounds.Qdotdots.upper'*ones(1,d*N));
opti.set_initial(A_col, guess.Qdotdots_col');

%% OCP: collocation equations
% Define CasADi variables for static parameters
tfk         = MX.sym('tfk');
% Define CasADi variables for states
ak          = MX.sym('ak',NMuscle);
aj          = MX.sym('akmesh',NMuscle,d);
akj         = [ak aj];
FTtildek    = MX.sym('FTtildek',NMuscle);
FTtildej    = MX.sym('FTtildej',NMuscle,d);
FTtildekj   = [FTtildek FTtildej];
Qsk         = MX.sym('Qsk',nq.all);
Qsj         = MX.sym('Qsj',nq.all,d);
Qskj        = [Qsk Qsj];
Qdotsk      = MX.sym('Qdotsk',nq.all);
Qdotsj      = MX.sym('Qdotsj',nq.all,d);
Qdotskj     = [Qdotsk Qdotsj];
a_ak        = MX.sym('a_ak',nq.arms);
a_aj        = MX.sym('a_akmesh',nq.arms,d);
a_akj       = [a_ak a_aj];
a_mtpk      = MX.sym('a_mtpk',nq.mtp);
a_mtpj      = MX.sym('a_mtpkmesh',nq.mtp,d);
a_mtpkj     = [a_mtpk a_mtpj];

% Define CasADi variables for controls
vAk     = MX.sym('vAk',NMuscle);
e_ak    = MX.sym('e_ak',nq.arms);
e_mtpk  = MX.sym('e_mtpk',nq.mtp);

% define the exoskeleton assistive torque
Texok   = MX.sym('Texo',2); % joint moments for the exoskeleton

% Define CasADi variables for "slack" controls
dFTtildej   = MX.sym('dFTtildej',NMuscle,d);
Aj          = MX.sym('Aj',nq.all,d);
eq_constr   = {}; % Initialize equality constraint vector
ineq_constr1 = {}; % Initialize inequality constraint vector 1
ineq_constr2 = {}; % Initialize inequality constraint vector 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time step
h = tfk/N;
% Loop over collocation points
for j=1:d
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Unscale variables
    Qskj_nsc = Qskj.*(scaling.QsQdots(1:2:end)'*ones(1,size(Qskj,2)/2));
    
    Qdotskj_nsc = Qdotskj.*(scaling.QsQdots(2:2:end)'* ...
        ones(1,size(Qdotskj,2)/2));
    FTtildekj_nsc = FTtildekj.*(scaling.FTtilde'*ones(1,size(FTtildekj,2)));
    dFTtildej_nsc = dFTtildej.*scaling.dFTtilde;
    Aj_nsc = Aj.*(scaling.Qdotdots'*ones(1,size(Aj,2)));
    vAk_nsc = vAk.*scaling.vA;
    
    QsQdotskj_nsc = MX(nq.all*2, d+1);
    QsQdotskj_nsc(1:2:end,:) = Qskj_nsc;
    QsQdotskj_nsc(2:2:end,:) = Qdotskj_nsc;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get muscle-tendon lengths, velocities, and moment arms
    % Left leg
    qinj_l          = Qskj_nsc(IndexLeft, j+1);
    qdotinj_l       = Qdotskj_nsc(IndexLeft, j+1);
    [lMTj_l,vMTj_l,MAj_l] =  f_lMT_vMT_dM(qinj_l,qdotinj_l);
    MAj.hip_flex.l   =  MAj_l(mai(1).mus.l',1);
    MAj.hip_add.l    =  MAj_l(mai(2).mus.l',2);
    MAj.hip_rot.l    =  MAj_l(mai(3).mus.l',3);
    MAj.knee.l       =  MAj_l(mai(4).mus.l',4);
    MAj.ankle.l      =  MAj_l(mai(5).mus.l',5);
    MAj.subt.l       =  MAj_l(mai(6).mus.l',6);
    % For the back muscles, we want left and right together: left
    % first, right second. In MuscleInfo, we first have the right
    % muscles (44:46) and then the left muscles (47:49). Since the back
    % muscles only depend on back dofs, we do not care if we extract
    % them "from the left or right leg" so here we just picked left.
    MAj.trunk_ext    =  MAj_l([47:49,mai(8).mus.l]',8);
    MAj.trunk_ben    =  MAj_l([47:49,mai(9).mus.l]',9);
    MAj.trunk_rot    =  MAj_l([47:49,mai(10).mus.l]',10);
    % Right leg
    qinj_r      = Qskj_nsc(IndexRight,j+1);
    qdotinj_r   = Qdotskj_nsc(IndexRight,j+1);
    [lMTj_r,vMTj_r,MAj_r] = f_lMT_vMT_dM(qinj_r,qdotinj_r);
    % Here we take the indices from left since the vector is 1:49
    MAj.ankle.r      =  MAj_r(mai(5).mus.l',5);
    MAj.subt.r       =  MAj_r(mai(6).mus.l',6);
    % Both legs
    % In MuscleInfo, we first have the right back muscles (44:46) and
    % then the left back muscles (47:49). Here we re-organize so that
    % we have first the left muscles and then the right muscles.
    lMTj_lr = [lMTj_l([1:43,47:49],1);lMTj_r(1:46,1)];
    vMTj_lr = [vMTj_l([1:43,47:49],1);vMTj_r(1:46,1)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get muscle-tendon forces and derive Hill-equilibrium
    [Hilldiffj,FTj,Fcej,Fpassj,Fisoj] = ...
        f_forceEquilibrium_FtildeState_all_tendon(akj(:,j+1),...
        FTtildekj_nsc(:,j+1),dFTtildej_nsc(:,j),...
        lMTj_lr,vMTj_lr,tensions);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get passive joint torques
    Tau_passj_all = f_AllPassiveTorques(Qskj_nsc(:,j+1),Qdotskj_nsc(:,j+1));
    Tau_passj.ankle.r = Tau_passj_all(10);
    Tau_passj.subt.r = Tau_passj_all(12);
    Tau_passj.tmt.r = Tau_passj_all(14);
    Tau_passj.mtp.r = Tau_passj_all(16);
    
    % Expression for the state derivatives at the collocation points
    Qsp_nsc      = Qskj_nsc*C(:,j+1);
    Qdotsp_nsc   = Qdotskj_nsc*C(:,j+1);
    FTtildep_nsc = FTtildekj_nsc*C(:,j+1);
    ap           = akj*C(:,j+1);
    % Append collocation equations
    % Dynamic constraints are scaled using the same scale
    % factors as the ones used to scale the states
    % Activation dynamics (implicit formulation)
    eq_constr{end+1} = (h*vAk_nsc - ap)./scaling.a;
    % Contraction dynamics (implicit formulation)
    eq_constr{end+1} = (h*dFTtildej_nsc(:,j) - FTtildep_nsc)./...
        scaling.FTtilde';
    % Skeleton dynamics (implicit formulation)
    qdotj_nsc = Qdotskj_nsc(:,j+1); % velocity
    eq_constr{end+1} = (h*qdotj_nsc - Qsp_nsc)./scaling.QsQdots(1:2:end)';
    eq_constr{end+1} = (h*Aj_nsc(:,j) - Qdotsp_nsc)./...
        scaling.QsQdots(2:2:end)';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Call external function (run inverse dynamics)
    [Tj] = F([QsQdotskj_nsc(:,j+1);Aj_nsc(:,j)]);    % left and right leg exoskeleton torques as inputs as well.

   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Knee, right
    eq_constr{end+1} = Tj(jointfi.tibia.ty,1) + F_tib_y; % vertical force on knee
    
    % Ankle, right
    Ft_ankle_r      = FTj(mai(5).mus.r',1);
    T_ankle_r       = f_T12(MAj.ankle.r,Ft_ankle_r);
    eq_constr{end+1} = Tj(jointfi.ankle.r,1)-(T_ankle_r + Tau_passj.ankle.r);
    
    % Subtalar, right
    Ft_subt_r       = FTj(mai(6).mus.r',1);
    T_subt_r        = f_T12(MAj.subt.r,Ft_subt_r);
    eq_constr{end+1} = Tj(jointfi.subt.r,1)-(T_subt_r + Tau_passj.subt.r );
    
    % Tmt, right
    eq_constr{end+1} = Tj(jointfi.tmt.r,1) - Tau_passj.tmt.r;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Activation dynamics (implicit formulation)
    act1 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tdeact);
    act2 = vAk_nsc + akj(:,j+1)./(ones(size(akj(:,j+1),1),1)*tact);
    ineq_constr1{end+1} = act1;
    ineq_constr2{end+1} = act2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Contraction dynamics (implicit formulation)
    eq_constr{end+1} = Hilldiffj;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
end % End loop over collocation points
eq_constr = vertcat(eq_constr{:});
ineq_constr1 = vertcat(ineq_constr1{:});
ineq_constr2 = vertcat(ineq_constr2{:});


% Casadi function to get constraints and objective
f_coll = Function('f_coll',{tfk,ak,aj,FTtildek,FTtildej,Qsk,Qsj,Qdotsk,...
    Qdotsj,a_ak,a_aj,a_mtpk,a_mtpj,vAk,e_ak,e_mtpk,dFTtildej,Aj,Texok},...
    {eq_constr,ineq_constr1,ineq_constr2});
% assign NLP problem to multiple cores
f_coll_map = f_coll.map(N,S.parallelMode,S.NThreads);
[coll_eq_constr, coll_ineq_constr1, coll_ineq_constr2] = f_coll_map(tf,...
    a(:,1:end-1), a_col, FTtilde(:,1:end-1), FTtilde_col, Qs(:,1:end-1), ...
    Qs_col, Qdots(:,1:end-1), Qdots_col, a_a(:,1:end-1), a_a_col, ...
    a_mtp(:,1:end-1), a_mtp_col, vA, e_a, e_mtp, dFTtilde_col, A_col,ExoVect);
% constrains
opti.subject_to(coll_eq_constr == 0);
opti.subject_to(coll_ineq_constr1(:) >= 0);
opti.subject_to(coll_ineq_constr2(:) <= 1/tact);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over mesh points
for k=1:N
    % Variables within current mesh interval
    % States
    akj = [a(:,k), a_col(:,(k-1)*d+1:k*d)];
    FTtildekj = [FTtilde(:,k), FTtilde_col(:,(k-1)*d+1:k*d)];
    Qskj = [Qs(:,k), Qs_col(:,(k-1)*d+1:k*d)];
    Qdotskj = [Qdots(:,k), Qdots_col(:,(k-1)*d+1:k*d)];
    a_akj = [a_a(:,k), a_a_col(:,(k-1)*d+1:k*d)];
    a_mtpkj = [a_mtp(:,k), a_mtp_col(:,(k-1)*d+1:k*d)];
    % Add equality constraints (next interval starts with end values of
    % states from previous interval)
    opti.subject_to(a(:,k+1) == akj*D);
    opti.subject_to(FTtilde(:,k+1) == FTtildekj*D); % scaled
    opti.subject_to(Qs(:,k+1) == Qskj*D); % scaled
    opti.subject_to(Qdots(:,k+1) == Qdotskj*D); % scaled
end % End loop over mesh points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional path constraints

    opti.subject_to(Qs(:,end) - Qs(:,1) == 0);
    opti.subject_to(Qdots(:,end) - Qdots(:,1) == 0);
    opti.subject_to(Qs(:,end) + Qs(:,1) == 0);
    opti.subject_to(Qdots(:,end) + Qdots(:,1) == 0);
    % Muscle activations
    opti.subject_to(a(:,end) - a(:,1) == 0);
    % Muscle-tendon forces
    opti.subject_to(FTtilde(:,end) - FTtilde(:,1) == 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create NLP solver
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.mu_strategy           = 'adaptive';
options.ipopt.max_iter              = S.max_iter;
options.ipopt.linear_solver         = S.linear_solver;
options.ipopt.tol                   = 1*10^(-S.tol_ipopt);
opti.solver('ipopt', options);
% Create and save diary
OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
if ~isfolder(OutFolder)
    mkdir(OutFolder);
end
Outname = fullfile(OutFolder,[S.savename '_log.txt']);
diary(Outname);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve problem
% Opti does not use bounds on variables but constraints. This function
% adjusts for that.
%     opti.solve();
[w_opt,stats] = solve_NLPSOL(opti,options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diary off
% Extract results
% Create setup
setup.tolerance.ipopt = S.tol_ipopt;
setup.bounds = bounds;
setup.scaling = scaling;
setup.guess = guess;

%% Save the results
Outname = fullfile(OutFolder,[S.savename '.mat']);
Sopt = S;
save(Outname,'w_opt','stats','setup','Sopt','ExoVect');
end

