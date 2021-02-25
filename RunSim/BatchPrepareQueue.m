clear all
close all
clc

%% Paths
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/OCP']);
addpath([pathRepo '/VariousFunctions']);
addpath([pathRepo '/Plots']);
addpath([pathRepo '/CasADiFunctions']);
addpath([pathRepo '/Musclemodel']);
addpath([pathRepo '/Polynomials']);

%% Manual settings
kTMT = [500 800 1000 1500]; % 250 2000
dTMT = [0 0.5]; % 0.2

exo = [[0; 0], [1; 0], [1; 1]]';

cWL = [0.02,0.03, 0.04];

fst=1;
for ik=1:length(kTMT)
    for id=1:length(dTMT)
        for ie=1:size(exo,1)
            for iw=1:length(cWL)
                
% tarsometatarsal joint
S.tmt = 1;              % 1: use a model with tmt joint
S.tmt_locked = 0;
% linear spring
S.kTMT = kTMT(ik);          % (Nm/rad) stiffness of tmt joint 
S.dTMT = dTMT(id);            % (Nms/rad) damping of tmt joint

% windlass mechanism
S.Windlass = 1;
S.cWL = cWL(iw);           % relative change in foot arch length at mtp 20° dorsiflexion

% assumption to simplify Hill-type muscle model
S.MuscModelAsmp = 0;    % 0: musc height = cst, 1: pennation angle = cst

% exo
S.ExoBool       = exo(ie,1);
S.ExoScale      = exo(ie,2);
                        
S.DataSet = 'PoggenSee2020_AFO';            % dataset with exoskeleton torque profile

% S.ExoImplementation = 'TorqueTibiaMetatarsi';
S.ExoImplementation = 'TorqueTibiaCalcn';

% Ideal assistance
ia = 0;
% S.T_max_ankle_exo = 30;
% S.P_max_ankle_exo = 50;

% output folder
S.ResultsFolder = 'batch_windlass';
suffixCasName = '';
suffixName = '';

% Folder with default functions
S.subject            = 's1_Poggensee';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 4;        % initial guess mode identifier (1 walk, 2 run, 3prev.solution, 4 solution from /IG/Data folder)

if S.IGmodeID == 4
S.savename_ig   = 'NoExo';
elseif S.IGmodeID == 3
S.ResultsF_ig   = 'PredSim_adaptations';
S.savename_ig   = '';
end

% Ideal assistance
if ia
    S.ExoController = 'Ideal Assistance';
end

% select folder with polynomials
S.PolyFolder = 's1_Poggensee';

% external function
if S.ExoBool == 0
    S.ExternalFunc  = 'PredSim_3D_Pog_s1_tmt_v3.dll';        % external function
    S.ExternalFunc2 = 'PredSim_3D_Pog_s1_tmt_pp_v3.dll';     % external function for post-processing
else
    S.ExternalFunc  = 'SimExo_3D_Pog_s1_tmt_TTC_v3.dll';
    S.ExternalFunc2  = 'SimExo_3D_Pog_s1_tmt_TTC_pp_v3.dll';
end

% build standardised names
[savename, casfuncfol] = getSavename(S);
S.CasadiFunc_Folders = [casfuncfol suffixCasName];
S.savename = [savename suffixName];

% if fst
    fst = 0;
    if (exist([pathRepo '/Results/batchQ.mat'],'file')==2) 
        load([pathRepo '/Results/batchQ.mat'],'batchQ');
    else
        batchQ.(S.savename) = struct('S',[]);
    end
% end
batchQ.(S.savename).S = S;

batchQ.(S.savename).PredSim = 'f_PredSim_Gait92_tmt';
batchQ.(S.savename).LoadSim = 'f_LoadSim_Gait92_tmt';

save([pathRepo '/Results/batchQ.mat'],'batchQ');
            end
        end
    end
end



