function [filteredResults]=getResultsForSameParams(dpath,varargin)

params = varargin{1,1};
nP = length(params);
% nP = length(varargin);
% params{1} = 'xxx';
% if nP>0
%     params{nP};
%     for i=1:nP
%         params{i} = varargin{i};
%     end
% end
% params{1} = param1;

% Get the names of the results files
MatFiles = dir(fullfile(dpath,'*.mat'));
nFil = length(MatFiles);
ct = 1;
for i = 1:nFil
    filename = MatFiles(i).name;
    OutName = fullfile(dpath,filename);
    if strcmp(filename(end-6:end),'_pp.mat')
        Names{ct,1} = OutName;
        Names{ct,2} = filename(1:end-7);
        ct= ct+1;
    end
end

nFil = ct-1;
idx = zeros(nP,nFil);
for i=1:nP
    for j=1:nFil
        if  length(params{i})>4 && strcmp(params{i}(1:4),'not_')
            bool = 0;
            param = params{i}(5:end);
        else
            bool = 1;
            param = params{i};
        end
        if sum(strfind(Names{j,2},param)) > 0
            idx(i,j)=bool;
        else
            idx(i,j)=abs(1-bool);
        end
    end
    if idx(i,:)==zeros(1,nFil)
        disp(['No files in the folder satisfy ' params{i} '. Dropping criterium.'])
        idx(i,:) = ones(1,nFil);
    end
end
ctr = 1;
filteredResults = {};
for j=1:nFil
    if sum(idx(:,j)) == nP
        filteredResults{ctr} = Names{j,1};
        ctr=ctr+1;
    end
end

end