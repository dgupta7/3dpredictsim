function [] = Plot3D_pwd(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dpath = pwd;
md = 1;
if ~isempty(varargin)
    dpath = varargin{1};
    if length(varargin)>1 && strcmp(varargin{2},'no_meas_data')
        md = 0;
    end
end


% Get the names of the results files
MatFiles = dir(fullfile(dpath,'*.mat'));
nFil = length(MatFiles);
ct = 1;
for i = 1:nFil
    filename = MatFiles(i).name;
    FileEnd = filename(end-6:end);
    OutName = fullfile(dpath,[filename(1:end-4) '_pp.mat']);
    if ~strcmp(FileEnd,'_pp.mat')
        Names{ct} = OutName;
        ct= ct+1;
    end
end

% plot influence Obj A
disp('Plotting: ');
DispHeader(Names);
disp('....');
CsV = hsv(ct);
h = figure();
set(h,'Position',[82         151        1497         827]);
for i = 1:ct-1
    [path,name,ext] = fileparts(Names{i});
    if md
        PlotResults_3DSim_tmt(Names{i},CsV(i,:),name,h);
    else
        PlotResults_3DSim_tmt(Names{i},CsV(i,:),name,h,varargin{2});
    end

end

% filename = fullfile(dpath,'FigureResults.fig');
% if ~exist(filename,'file')
% saveas(h,fullfile(dpath,'FigureResults.fig'));
% end


end

