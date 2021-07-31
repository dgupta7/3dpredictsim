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
addpath([pathRepo '/Debug']);
addpath([pathRepo '/FootModel']);
AddCasadiPaths();

%% Manual settings

% kMT = [30, 50, 100, 150, 200, 250, 300, 400, 500, 800, 1000, 1500, 2000];
% kMT = [50, 100, 150, 200, 250, 300, 400, 500, 800];
% kMT = [250, 500, 800, 1000, 1500, 2000];
kMT = [5000, 1e4, 1e5, 1e6];
dMT = [0]; % 0,1,2,5
v = 1.25;

% kMT = 300;
% dMT = 0;
% v = [0.8, 1, 1.2, 1.4, 1.6, 1.8, 2.7];



fst=1;
for ik=1:length(kMT)
    for id=1:length(dMT)
        for iv=1:length(v)
                
% no tmt
S.tmt = 0; 

% with mtj
S.mtj = 1; 
% plantar fascia
S.PF_stiffness = 'Natali2010'; %'Natali2010'
S.sf_PF = 1;
S.PF_slack_length = 0.15;
S.R_mtth = 9.5e-3;
S.WLpoly = 1;
% mtj
S.MT_li_nonl = 0;
S.mtj_stiffness = 'Song2011';
S.kMT_li = kMT(ik);
S.kMT_li2 = 10;
S.dMT = dMT(id);
% mtp
S.WL_T_mtp = 1;
S.Mu_mtp = 0; 
S.kMTP = 5;
% S.dMTP = 0;

% assumption to simplify Hill-type muscle model
S.MuscModelAsmp = 0;    % 0: musc height = cst, 1: pennation angle = cst
S.contactStiff = 1;
S.v_tgt     = v(iv);     % average speed 1.25
S.N         = 50;       % number of mesh intervals

% exo
S.ExoBool       = 0;
S.ExoScale      = 0;          
S.DataSet = 'PoggenSee2020_AFO';            % dataset with exoskeleton torque profile

% S.ExoImplementation = 'TorqueTibiaMetatarsi';
S.ExoImplementation = 'TorqueTibiaCalcn';


% output folder
S.ResultsFolder = 'MidTarsalJoint';
suffixCasName = '';
% suffixName = ['_v' num2str(v(iv)*10)];
suffixName = '';

% Folder with default functions
S.subject = 'subject1';

% initial guess based on simulations without exoskeletons
S.IGsel         = 2;        % initial guess identifier (1: quasi random, 2: data-based)
S.IGmodeID      = 3;        % initial guess mode identifier (1 walk, 2 run, 3 prev.solution, 4 solution from /IG/Data folder)

if S.IGmodeID == 4
    S.savename_ig   = 'NoExo';
elseif S.IGmodeID == 3
    S.ResultsF_ig   = 'MuscleModel';
    if strcmp(S.subject,'s1_Poggensee')
        S.savename_ig   = 'Pog_s1_bCst_ig24';
    else
        S.savename_ig   = 'Fal_s1_bCst_ig1';
    end
end
if S.IGmodeID == 1
    if strcmp(S.subject,'s1_Poggensee')
        S.IG_PelvisY = 0.896 + 0.0131;
    else
        S.IG_PelvisY = 0.9385 + 0.0131;
    end
end

% select folder with polynomials
S.PolyFolder = S.subject;

% external function
% external function
if S.tmt == 0 && S.mtj == 0
    if strcmp(S.subject,'s1_Poggensee')
        if S.ExoBool == 0
            S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtp.dll';
            S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtp_pp.dll';
        else
            S.ExternalFunc  = 'SimExo_3D_talus_out.dll';
        end
    elseif strcmp(S.subject,'subject1')
        if S.ExoBool == 0
%             S.ExternalFunc  = 'ID_Subject1.dll';
%             S.ExternalFunc2 = 'PredSim_3D_Fal_s1_pp_v2.dll';
            S.ExternalFunc  = 'PredSim_3D_Fal_s1_v7.dll';
            S.ExternalFunc2 = 'PredSim_3D_Fal_s1_pp_v4.dll';
        end
    end
    
elseif S.mtj == 1
    if S.ExoBool == 0
        if strcmp(S.subject,'s1_Poggensee')
            S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtj_v3.dll';
            S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtj_pp_v3.dll';
            
        elseif strcmp(S.subject,'subject1')
             if S.contactStiff == 10
                S.ExternalFunc  = 'PredSim_3D_Fal_s1_mtj_spx10_v1.dll';
                S.ExternalFunc2 = 'PredSim_3D_Fal_s1_mtj_spx10_pp_v1.dll';
             elseif S.contactStiff == 2
                S.ExternalFunc  = 'PredSim_3D_Fal_s1_mtj_spx2_v1.dll';
                S.ExternalFunc2 = 'PredSim_3D_Fal_s1_mtj_spx2_pp_v1.dll';
             else
                S.ExternalFunc  = 'PredSim_3D_Fal_s1_mtj_v1.dll';
                S.ExternalFunc2 = 'PredSim_3D_Fal_s1_mtj_pp_v5.dll';
             end
        end
    end
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

if S.tmt == 0 && S.mtj == 0
    batchQ.(S.savename).PredSim = 'f_PredSim_Gait92';
    batchQ.(S.savename).LoadSim = 'f_LoadSim_Gait92';

else
    if S.Mu_mtp
        batchQ.(S.savename).PredSim = 'f_PredSim_Gait92_tmt_v2';
        batchQ.(S.savename).LoadSim = 'f_LoadSim_Gait92_tmt_v2';
    else
        batchQ.(S.savename).PredSim = 'f_PredSim_Gait92_tmt';
        batchQ.(S.savename).LoadSim = 'f_LoadSim_Gait92_tmt';
    end
end


save([pathRepo '/Results/batchQ.mat'],'batchQ');


        end
    end
end



