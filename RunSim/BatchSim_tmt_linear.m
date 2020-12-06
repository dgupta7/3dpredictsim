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

% settings for optimization
S.v_tgt     = 1.25;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 2;        % number of threads for parallel computing
S.max_iter  = 10000;    % maximum number of iterations

% tarsometatarsal joint
S.tmt = 1;              % 1: use a model with tmt joint
S.tmt_locked = 0;

kTMT = [250 500 800 1000 2000];
dTMT = [0 0.2 0.5];
exo = [[0; 0], [1; 0], [1; 1]]';


% set up settings fo
count = 1;
for ik=1:length(kTMT)
    for id=1:length(dTMT)
        for ie=1:3
            

S.kTMT = kTMT(ik);
S.dTMT = dTMT(id);

% assumption to simplify Hill-type muscle model
S.MuscModelAsmp = 0;    % 0: musc width = cst, 1: pennation angle = cst

% exo
S.ExoBool       = exo(ie,1);    % 1: is wearing exo
S.ExoScale      = exo(ie,2);    % scale factor of exoskeleton assistance profile 
                        % 0: no assistance (passive) 1: nominal assistance (active)
                        
S.DataSet = 'PoggenSee2020_AFO';            % dataset with exoskeleton torque profile
                        
% output folder
S.ResultsFolder = 'Batchsim_tmt_linear';
% S.ResultsFolder = 'debug_batch';

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


%% Automated settings
pathResults = fullfile([pathRepo '/Results'],S.ResultsFolder);
if ~isfolder(pathResults)
    mkdir(pathResults);
end

% build savename and casadifunction foldername
if strcmp(S.subject,'s1_Poggensee')
    savename = 'Pog_s1';
    casfuncfol = 'casadi_s1Pog';
end
if S.tmt
    savename = [savename '_tmt'];
    casfuncfol = [casfuncfol '_tmt'];
    if S.tmt_locked
        savename = [savename 'L'];
    end
end
if isfield(S,'MuscModelAsmp') && ~isempty(S.MuscModelAsmp) && S.MuscModelAsmp==0
    savename = [savename '_bCst'];
    casfuncfol = [casfuncfol '_MuscModel_bCst'];
else
    savename = [savename '_aCst'];
    casfuncfol = [casfuncfol '_MuscModel_alphaCst'];
end
if S.tmt && isfield(S,'dTMT') && ~isempty(S.dTMT) && ~S.tmt_locked
    savename = [savename '_d0' num2str(S.dTMT*10)];
    casfuncfol = [casfuncfol '_d0' num2str(S.dTMT*10)];
end
if S.tmt && isfield(S,'kTMT') && ~isempty(S.kTMT) && ~S.tmt_locked
    savename = [savename '_k' num2str(S.kTMT)];
    casfuncfol = [casfuncfol '_k' num2str(S.kTMT)];
end
if S.IGsel == 1
    savename = [savename '_ig1'];
else
    savename = [savename '_ig2' num2str(S.IGmodeID )];
end
if S.ExoBool == 1
    if S.ExoScale == 0
        savename = [savename '_pas'];
    else
        savename = [savename '_act'];
    end    
end

S.CasadiFunc_Folders = casfuncfol;
S.savename = savename;

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% select folder with polynomials
S.PolyFolder = 's1_Poggensee';

% external function
if S.ExoBool == 0
    S.ExternalFunc  = 'PredSim_3D_Pog_s1_tmt.dll';        % external function
    S.ExternalFunc2 = 'PredSim_3D_Pog_s1_tmt_pp.dll';     % external function for post-processing
else
    S.ExternalFunc  = 'SimExo_3D_Pog_s1_tmt.dll';
    S.ExternalFunc2  = 'SimExo_3D_Pog_s1_tmt_pp.dll';
end
    


% Create the casadifunctions if they do not exist yet
if ~isfolder([pathRepo '\CasADiFunctions\' S.CasadiFunc_Folders])
    CreateCasADiFunctions_all_tmt(pathRepo,S);
end

S_batch{count} = S;
count = count+1;
clear('savename','casfuncfol');
        end
    end
end

%% Run
name = getenv('COMPUTERNAME');
if strcmp(name,'GBW-D-W2711')
myCluster = parcluster('LocalProfile1_Lars_10x2');
elseif strcmp(name,'MSI')   %Lars
myCluster = parcluster('LocalProfile_2x2');
end

StartPath = pwd;
MainPath = pathRepo;

PathPolynomials = fullfile(MainPath,'Polynomials',S.subject);
ExoPath = fullfile(MainPath,'Data','Poggensee_2020');
pathExternalFunctions = fullfile(MainPath,'ExternalFunctions');
pathResults = fullfile([pathRepo '/Results'],S.ResultsFolder);

j=1;
imax = length(S_batch);
for i=1:imax
CasadiFiles = fullfile(MainPath,'CasADiFunctions',S_batch{i}.CasadiFunc_Folders);
pathResult_pp = fullfile([pathRepo '/Results'],S_batch{i}.ResultsFolder,[S_batch{i}.savename '_pp.mat']);
    if ~exist(pathResult_pp,'file')
        % solve
        job(j) = batch(myCluster,'f_PredSim_Gait92_tmt',0,{S_batch{i}},'CurrentFolder',StartPath,'AdditionalPaths',{CasadiFiles,PathPolynomials,ExoPath,pathExternalFunctions});
        j=j+1;
        % % post-proces
        % job(j) = batch(myCluster,'f_LoadSim_Gait92_tmt',0,{pathResult,S_batch{i}.savename},'CurrentFolder',StartPath,'AdditionalPaths',{CasadiFiles,PathPolynomials,ExoPath,pathExternalFunctions});
        % j=j+1;
    end
end
%% rerun this section after the jobs are done to get logfiles

for i=1:length(job)
    if strcmp(job(1, i).State,'finished') && strcmp(job(1, i).Name,'f_PredSim_Gait92_tmt')
        diary(fullfile(pathRepo,'Results',S_batch{i}.ResultsFolder,[S_batch{i}.savename '_log.txt']));
        job(1, i).Tasks.Diary
        diary off
    end
end
