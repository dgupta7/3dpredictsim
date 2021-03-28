function [filteredResults]=filterResultfolderByParameters(dpaths,params)
% Folder will be filtered to only return results that satisfy all chosen
% parameter settings. If a parameter matches none of the results, it is not
% used as a filter. (To prevent wrong entries or typos.)


ctr = 1;
filteredResults = {};
    
for k=1:numel(dpaths)
    dpath = dpaths{k};
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
        elseif strcmp(filename(1:5),'Foot_')
            Names{ct,1} = OutName;
            Names{ct,2} = filename(1:end-4);
            ct= ct+1;
        end
    end

    nP = length(params);
    nN = ct-1;
    idx = zeros(nP,nN);
    for i=1:nP
        for j=1:nN
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
        if idx(i,:)==zeros(1,nN)
            disp(['No files in the folder "' dpath '" satisfy ' params{i} ', dropping criterium.'])
            idx(i,:) = ones(1,nN);
        end
    end
    
    for j=1:nN
        if sum(idx(:,j)) == nP
            filteredResults{ctr} = Names{j,1};
            ctr=ctr+1;
        end
    end

end

end