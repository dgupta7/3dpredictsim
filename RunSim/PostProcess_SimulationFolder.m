%% PostProcess Simluations
%---------------------------


% Default script for post-processing of all simulations
%------------------------------------------------------
clear all; clc;


%% Path information
cd ..;
Datapath = [pwd '\Results'];
addpath([pwd,'/OCP']);
addpath([pwd,'/MuscleModel']);
addpath([pwd,'/Debug']);
addpath([pwd '/VariousFunctions']);
AddCasadiPaths();
DataFolders = {'windlass'};


S.OverWrite = 0;

%% Post process the simulations

ct = 1;
nF = length(DataFolders);
for f = 1:nF
    % get all the simulation results in this folder
    dpath = fullfile(Datapath,DataFolders{f});
    MatFiles = dir(fullfile(dpath,'*.mat'));
    nFil = length(MatFiles);
    for i = 1:nFil
        filename = MatFiles(i).name;
        FileEnd = filename(end-6:end);
        OutName = fullfile(dpath,[filename(1:end-4) '_pp.mat']);
        if ~strcmp(FileEnd,'_pp.mat')
            Names{ct} = OutName;
            FolderIndex(ct) = f;
            ct= ct+1;
            if (~exist(OutName,'file') || S.OverWrite == 1)
                try
                    f_LoadSim_Gait92_tmt(DataFolders{f},filename);
                catch
                    f_LoadSim_Gait92(DataFolders{f},filename);
                end
%                 disp(ct);
            end
        end
    end
end