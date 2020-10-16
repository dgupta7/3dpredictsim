[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);

dpath = [pathRepo,'/Results'];
OutPath = [pathRepo,'/Figures'];
addpath([pathRepo '/VariousFunctions']);

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