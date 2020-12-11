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
S.tmt = 0;              % 1: use a model with tmt joint
S.tmt_locked = 0;

% assumption to simplify Hill-type muscle model
mumo = [0 1];    % 0: musc height = cst, 1: pennation angle = cst

exo = [[0; 0], [1; 0], [1; 1]]';


count = 1;
for ie=1:size(exo,1)
    for im=1:length(mumo)

% assumption to simplify Hill-type muscle model
S.MuscModelAsmp = mumo(im);

% exo
S.ExoBool       = exo(ie,1);    % 1: is wearing exo
S.ExoScale      = exo(ie,2);    % scale factor of exoskeleton assistance profile 
                        % 0: no assistance (passive) 1: nominal assistance (active)
                        
S.DataSet = 'PoggenSee2020_AFO';            % dataset with exoskeleton torque profile
                        
% output folder
S.ResultsFolder = 'Batchsim_without_tmt';
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

[savename, casfuncfol] = getSavename(S);
S.CasadiFunc_Folders = casfuncfol;
S.savename = savename;

% quasi random initial guess, pelvis y position
S.IG_PelvisY = 0.896;   % subject 1 poggensee

% select folder with polynomials
S.PolyFolder = 's1_Poggensee';

% external function
if S.ExoBool == 0
    S.ExternalFunc  = 'PredSim_3D_Pog_s1_mtp.dll';        % external function
    S.ExternalFunc2 = 'PredSim_3D_Pog_s1_mtp_pp.dll';     % external function for post-processing
else
    S.ExternalFunc  = 'SimExo_3D_talus_out.dll';        % this one is with the pinjoint mtp
end
    


% Create the casadifunctions if they do not exist yet
if ~isfolder([pathRepo '\CasADiFunctions\' S.CasadiFunc_Folders])
    disp('Creating casadifunctions...');
    CreateCasADiFunctions_all_tmt(pathRepo,S);
    disp('...casadifunctions created');
end

S_batch{count} = S;
count = count+1;
clear('savename','casfuncfol');

    end
end


%% Run
name = getenv('COMPUTERNAME');
if strcmp(name,'GBW-D-W2711')
myCluster = parcluster('LocalProfile1_Lars_8x2');
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
            job(j) = batch(myCluster,'f_PredSim_Gait92',0,{S_batch{i}},'CurrentFolder',StartPath,...
                'AdditionalPaths',{CasadiFiles,PathPolynomials,ExoPath,pathExternalFunctions});
        j=j+1;
    end
end

%% rerun this section after the jobs are done to get logfiles

for i=1:length(job)
    if strcmp(job(1, i).State,'finished') && strcmp(job(1, i).Name(1:9),'f_PredSim')
        diary(fullfile(pathRepo,'Results',S_batch{i}.ResultsFolder,[S_batch{i}.savename '_log.txt']));
        job(1, i).Tasks.Diary
        diary off
    end
end

