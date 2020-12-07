function [] = Plot3D_pwd_separate(varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dpath = pwd;
md = 1;
if ~isempty(varargin)
    dpath = varargin{1};
    if length(varargin)>1 
        if strcmp(varargin{2},'no_meas_data')
            md = 0;
            meas_type = {'no_meas_data'};
        end
    end
end
pathData = 'D:\school\WTK\thesis\model\3dpredictsim/ExperimentalData/ExperimentalData.mat';

% Get the names of the results files
MatFiles = dir(fullfile(dpath,'*.mat'));
nFil = length(MatFiles);
ct = [1,1,1];
for i = 1:nFil
    filename = MatFiles(i).name;
    FileEnd = filename(end-6:end);
    OutName = fullfile(dpath,[filename(1:end-4) '_pp.mat']);
    if ~strcmp(FileEnd,'_pp.mat')
        if strcmp(FileEnd,'act.mat')
            Names_act{ct(1)} = OutName;
            ct(1)= ct(1)+1;
        elseif strcmp(FileEnd,'pas.mat')
            Names_pas{ct(2)} = OutName;
            ct(2)= ct(2)+1;
        else
            Names{ct(3)} = OutName;
            ct(3)= ct(3)+1;
        end
    end
end

% plot w/o exo
if md
    meas_type = {'norm'};
end
disp('Plotting: ');
DispHeader(Names);
disp('....');
CsV = hsv(ct(3));
h = figure();
set(h,'Position',[82         151        1497         827]);
for i = 1:ct(3)-1
    [path,name,ext] = fileparts(Names{i});
    PlotResults_3DSim_tmt(Names{i},CsV(i,:),name,h,meas_type);
end
ValidationPlots(pathData,Names{1},Names{2:end})

% plot passive exo
if md
    meas_type = {'pas'};
end
disp('Plotting: ');
DispHeader(Names);
disp('....');
CsV = hsv(ct(2));
h_pas = figure();
set(h_pas,'Position',[82         151        1497         827]);
for i = 1:ct(2)-1
    [path,name,ext] = fileparts(Names_pas{i});
    PlotResults_3DSim_tmt(Names_pas{i},CsV(i,:),name,h_pas,meas_type);
end
ValidationPlots(pathData,Names_pas{1},Names_pas{2:end})

% plot active exo
if md
    meas_type = {'act'};
end
disp('Plotting: ');
DispHeader(Names);
disp('....');
CsV = hsv(ct(1));
h_act = figure();
set(h_act,'Position',[82         151        1497         827]);
for i = 1:ct(1)-1
    [path,name,ext] = fileparts(Names_act{i});
    PlotResults_3DSim_tmt(Names_act{i},CsV(i,:),name,h_act,meas_type);
end
ValidationPlots(pathData,Names_act{1},Names_act{2:end})

% filename = fullfile(dpath,'FigureResults.fig');
% if ~exist(filename,'file')
% saveas(h,fullfile(dpath,'FigureResults.fig'));
% end


end

