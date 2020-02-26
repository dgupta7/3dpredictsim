
%% Batch Sim exoskeleton assistance

% Default simulations antoine with exoskeleton assistance.
% Note currently without inertia of the exoskeleton
clear all;

%% Default settings

% flow control
S.Flow.solveProblem     = 1;   % set to 1 to solve problem
S.Flow.analyseResults   = 1;   % set to 1 to analyze results
S.Flow.loadResults      = 0;   % set to 1 to load results
S.Flow.saveResults      = 1;   % set to 1 to save sens. results
S.Flow.checkBoundsIG    = 0;   % set to 1 to visualize guess-bounds
S.Flow.writeIKmotion    = 1;   % set to 1 to write .mot file

% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.W.E       = 500;      % weight metabolic energy rate
S.W.Ak      = 50000;    % weight joint accelerations
S.W.ArmE    = 1000000;  % weight arm excitations
S.W.passMom = 1000;     % weight passive torques
S.W.A       = 2000;     % weight muscle activations
S.W.exp_E   = 2;        % power metabolic energy
S.W.Mtp     = 1000;     % weight mtp excitations
S.W.u       = 0.001;    % weight on excitations arms actuators
S.IGsel     = 2;        % initial guess identifier
S.ContactID = 1;        % contact model identifier
S.IGmodeID  = 1;        % initial guess mode identifier
S.IGcase    = 0;        % initial guess case identifier
S.h_weak    = 0;        % weakness hip actuators
S.Max_s     = 0;        % maximal contraction velocity identifier
S.pf_weak   = 0;        % weakness ankle plantaflexors
S.mE        = 0;        % metabolic energy model identifier
S.coCont    = 0;        % co-contraction identifier

% save name initial guess
S.savename_ig =[];

% ipopt options
S.linear_solver = 'mumps';
S.tol_ipopt     = 4;

% external function 
S.ExternalFunc = 'PredSim_mtp_cm1.dll';
S.ExternalFunc2 = 'PredSim_mtp_pp_cm1.dll';

% Folder with default functions
S.CasadiFunc_Folders = 'Casadi_Default';


%% Specific settings to runs simulation

% Simulation without exoskeleton assistance
% S.savename      = 'NoAssistance';
% S.loadname      = 'NoAssistance';
% S.ResultsFolder = 'IdealAssistance_NoInertia';
% S.ExoBool       = 0;
% S.ExoScale      = 0;
% S.DataSet       = 'PoggenSee2020_AFO';
% RP = f_PredSim_PoggenSee2020_DefaultS(S);

% Simulation with 50% assistance
S.savename      = 'Assistance50';
S.loadname      = 'Assistance50';
S.ResultsFolder = 'IdealAssistance_NoInertia';
S.ExoBool       = 1;
S.ExoScale      = 0.5;
S.DataSet       = 'PoggenSee2020_AFO';
R1 = f_PredSim_PoggenSee2020_DefaultS(S);

% Simulation with 100% assistance
S.savename      = 'Assistance100';
S.loadname      = 'Assistance100';
S.ResultsFolder = 'IdealAssistance_NoInertia';
S.ExoBool       = 1;
S.ExoScale      = 1;
S.DataSet       = 'PoggenSee2020_AFO';
R2 = f_PredSim_PoggenSee2020_DefaultS(S);

% Simulation with 150% assistance
S.savename      = 'Assistance150';
S.loadname      = 'Assistance150';
S.ResultsFolder = 'IdealAssistance_NoInertia';
S.ExoBool       = 1;
S.ExoScale      = 1.5;
S.DataSet       = 'PoggenSee2020_AFO';
R3 = f_PredSim_PoggenSee2020_DefaultS(S);

% Simulation with 200% assistance
S.savename      = 'Assistance200';
S.loadname      = 'Assistance200';
S.ResultsFolder = 'IdealAssistance_NoInertia';
S.ExoBool       = 1;
S.ExoScale      = 2;
S.DataSet       = 'PoggenSee2020_AFO';
R4 = f_PredSim_PoggenSee2020_DefaultS(S);

% Simulation with 250% assistance
S.savename      = 'Assistance250';
S.loadname      = 'Assistance250';
S.ResultsFolder = 'IdealAssistance_NoInertia';
S.ExoBool       = 1;
S.ExoScale      = 2.5;
S.DataSet       = 'PoggenSee2020_AFO';
R4 = f_PredSim_PoggenSee2020_DefaultS(S);


% Simulation with 300% assistance
S.savename      = 'Assistance300';
S.loadname      = 'Assistance300';
S.ResultsFolder = 'IdealAssistance_NoInertia';
S.ExoBool       = 1;
S.ExoScale      = 3;
S.DataSet       = 'PoggenSee2020_AFO';
R4 = f_PredSim_PoggenSee2020_DefaultS(S);