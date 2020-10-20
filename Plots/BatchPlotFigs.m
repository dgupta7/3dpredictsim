clear
clc

[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);

dpath = [pathRepo,'/Results'];
OutPath = [pathRepo,'/Figures'];
addpath([pathRepo '/VariousFunctions']);

%% plot batch
FolderNames = {'Test_Lars'};
OutFNames = {'Test_Lars'};

nf= length(FolderNames);

for i= 1:nf
    Fsel = fullfile(dpath,FolderNames{i});
    Fout = fullfile(OutPath,OutFNames{i});
    Plot3D_pwd(Fsel);
    OutFile = fullfile(Fout,'FigureResults.fig');
    if exist(OutFile,'file')
        delete(OutFile);
    end
    copyfile(fullfile(Fsel,'FigureResults.fig'),fullfile(Fout,'FigureResults.fig'));
end

%% plot difference
% addpath([dpath,'/Test_Lars']);
% ResultsFile1 = fullfile(dpath,'Test_Lars','Test1_bCst_exo_pp.mat');
% ResultsFile2 = fullfile(dpath,'Test_Lars','Test1_alphaCst_exo_pp.mat');
% PlotResultsComparison_3DSim(ResultsFile1,ResultsFile2)
