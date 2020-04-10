function [] = PostProcess_pwd(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dpath = pwd;
if ~isempty(varargin)
    dpath = varargin{1};
end


% Get the names of the results files
MatFiles = dir(fullfile(dpath,'*.mat'));
nFil = length(MatFiles);
ct = 1;
for i = 1:nFil
    filename = MatFiles(i).name;
    FileEnd = filename(end-6:end);
    OutName = fullfile(dpath,filename);
    if ~strcmp(FileEnd,'_pp.mat')
        Names{ct} = OutName;
        ct= ct+1;
    end
end

% plot influence Obj A
disp('PostProcessing: ');
DispHeader(Names);
disp('....');
CsV = hsv(ct);
for i = 1:ct-1
    [path,name,ext] = fileparts(Names{i});
    [mainR,Fname ]= fileparts(path);
    f_LoadSim_PoggenSee2020_DefaultS(Fname,name);
end



end

