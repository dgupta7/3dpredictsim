function [savename, casfuncfol,varargout] = getSavename(S)
% build savename and casadifunction foldername
savenameparts = {};
casfuncfolparts = {};
not = {};

% general settings
savenameparts{end+1} = S.ExternalFunc;
casfuncfolparts{end+1} = S.ExternalFunc;

if S.AchillesTendonScaleFactor~=1
    savenameparts{end+1} = ['ATx' num2str(S.AchillesTendonScaleFactor*100)];
    casfuncfolparts{end+1} = ['ATx' num2str(S.AchillesTendonScaleFactor*100)];
end

% mtp related settings
savenameparts{end+1} = 'MTP';
casfuncfolparts{end+1} = 'MTP';
if S.Foot.mtp_muscles
    savenameparts{end} = [savenameparts{end} 'm'];
    casfuncfolparts{end} = [casfuncfolparts{end} 'm'];
elseif S.Foot.mtp_actuator
    savenameparts{end} = [savenameparts{end} 'a'];
    casfuncfolparts{end} = [casfuncfolparts{end} 'a'];
else
    savenameparts{end} = [savenameparts{end} 'p'];
    casfuncfolparts{end} = [casfuncfolparts{end} 'p'];
end

if isfield(S.Foot,'kMTP') && ~isempty(S.Foot.kMTP)
    savenameparts{end+1} = ['k' num2str(S.Foot.kMTP)];
    casfuncfolparts{end+1} = ['k' num2str(S.Foot.kMTP)];
end
if isfield(S.Foot,'dMTP') && ~isempty(S.Foot.dMTP)
    savenameparts{end+1} = ['d0' num2str(S.Foot.dMTP*10)];
    casfuncfolparts{end+1} = ['d0' num2str(S.Foot.dMTP*10)];
end


% mtj related settings
if strcmp(S.Foot.Model,'mtj')
    if S.Foot.mtj_muscles
        savenameparts{end+1} = ['MTJm'];
        casfuncfolparts{end+1} = ['MTJm'];
    else
        savenameparts{end+1} = ['MTJp'];
        casfuncfolparts{end+1} = ['MTJp'];
    end
    if isfield(S.Foot,'MT_li_nonl') && ~isempty(S.Foot.MT_li_nonl) && S.Foot.MT_li_nonl
        if isfield(S.Foot,'mtj_stiffness') && ~isempty(S.Foot.mtj_stiffness)
            if strcmp(S.Foot.mtj_stiffness,'signed_lin')
                savenameparts{end+1} = ['nl_k' num2str(S.Foot.kMT_li) '_' num2str(S.Foot.kMT_li2)];
                casfuncfolparts{end+1} = ['nl_k' num2str(S.Foot.kMT_li) '_' num2str(S.Foot.kMT_li2)];
            else
                savenameparts{end+1} = ['nl_' S.Foot.mtj_stiffness];
                casfuncfolparts{end+1} = ['nl_' S.Foot.mtj_stiffness];
            end
        else
            savenameparts{end+1} = ['nl'];
            casfuncfolparts{end+1} = ['nl'];
        end
    elseif isfield(S.Foot,'kMT_li') && ~isempty(S.Foot.kMT_li)
        savenameparts{end+1} = ['k' num2str(S.Foot.kMT_li)];
        casfuncfolparts{end+1} = ['k' num2str(S.Foot.kMT_li)];
    end
    if isfield(S.Foot,'MT_li_nonl') && ~isempty(S.Foot.MT_li_nonl) && ~S.Foot.MT_li_nonl
        not{end+1} = 'not_MTJp_nl';
        not{end+1} = 'not_MTJm_nl';
    end
    if isfield(S.Foot,'dMT') && ~isempty(S.Foot.dMT) && S.Foot.dMT~=0
        savenameparts{end+1} = ['d0' num2str(S.Foot.dMT*10)];
        casfuncfolparts{end+1} = ['d0' num2str(S.Foot.dMT*10)];
    end

    if isfield(S.Foot,'PF_stiffness') && ~isempty(S.Foot.PF_stiffness)
        savenameparts{end+1} = ['PF_' S.Foot.PF_stiffness];
        casfuncfolparts{end+1} = ['PF_' S.Foot.PF_stiffness];
    end
    if S.Foot.PF_sf~=1
        savenameparts{end} = [savenameparts{end} '_x' num2str(S.Foot.PF_sf)];
        casfuncfolparts{end} = [casfuncfolparts{end} '_x' num2str(S.Foot.PF_sf)];
    end
    if isfield(S.Foot,'PF_slack_length') && ~isempty(S.Foot.PF_slack_length)
        savenameparts{end+1} = ['ls' num2str(S.Foot.PF_slack_length*1000)];
        casfuncfolparts{end+1} = ['ls' num2str(S.Foot.PF_slack_length*1000)];
    end

    if S.Foot.PIM==1
        savenameparts{end+1} = 'PIM';
        casfuncfolparts{end+1} = 'PIM';
        if isfield(S.W,'PIM') && ~isempty(S.W.PIM)
            savenameparts{end+1} = ['w' num2str(S.W.PIM,2)];
        end
        if isfield(S.W,'P_PIM') && ~isempty(S.W.P_PIM)
            savenameparts{end+1} = ['w' num2str(S.W.P_PIM,2)];
        end
    end
end
    


if isfield(S,'IGsel') && ~isempty(S.IGsel)
    if S.IGsel == 1
        savenameparts{end+1} = ['ig1'];
    else
        savenameparts{end+1} = ['ig2' num2str(S.IGmodeID )];
    end
end

if isfield(S,'suffixCasName')
    casfuncfolparts{end+1} = S.suffixCasName;
end
if isfield(S,'suffixName')
    savenameparts{end+1} = S.suffixName;
end
savename = savenameparts{1};
for i=2:numel(savenameparts)
    savename = [savename '_' savenameparts{i}];
end
criteria = {};
for i=1:numel(savenameparts)
    criteria{end+1} = savenameparts{i};
end
casfuncfol = casfuncfolparts{1};
for i=2:numel(casfuncfolparts)
    casfuncfol = [casfuncfol '_' casfuncfolparts{i}];
end

for i=1:numel(not)
    criteria{end+1} = not{i};
end


if nargout == 3
    varargout{1} = criteria;
else
    disp(savename);
    disp(casfuncfol);
end

end