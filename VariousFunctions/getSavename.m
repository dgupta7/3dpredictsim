function [savename, casfuncfol,varargout] = getSavename(S)
% build savename and casadifunction foldername
savenameparts = {};
casfuncfolparts = {};
not = {};

if isfield(S,'subject') && ~isempty(S.subject)
    if strcmp(S.subject,'s1_Poggensee')
        savenameparts{end+1} = 'Pog_s1';
        casfuncfolparts{end+1} = 'casadi_s1Pog';
    elseif strcmp(S.subject,'subject1')
        savenameparts{end+1} = 'Fal_s1';
        casfuncfolparts{end+1} = 'casadi_s1Fal';
    end
end
if isfield(S,'tmt') && ~isempty(S.tmt)
    if S.tmt
        savenameparts{end+1} = ['tmt'];
        casfuncfolparts{end+1} = ['tmt'];
%         if S.tmt_locked
%             savenameparts{end} = [savenameparts{end} 'L'];
%             casfuncfolparts{end} = [casfuncfolparts{end} 'L'];
%         else
%             not{end+1} ='not_tmtL';
%         end
    else
        not{end+1} ='not_tmt';
    end
end
if isfield(S,'MuscModelAsmp') && ~isempty(S.MuscModelAsmp) && S.MuscModelAsmp==0
    savenameparts{end+1} = ['bCst'];
    casfuncfolparts{end+1} = ['MuscModel_bCst'];
else
    savenameparts{end+1} = ['aCst'];
    casfuncfolparts{end+1} = ['MuscModel_alphaCst'];
end
nr = numel(savenameparts);
if isfield(S,'tmt') && ~isempty(S.tmt)
    if S.tmt && isfield(S,'dTMT') && ~isempty(S.dTMT) && ~S.tmt_locked
        savenameparts{end+1} = ['d0' num2str(S.dTMT*10)];
        casfuncfolparts{end+1} = ['d0' num2str(S.dTMT*10)];
    end
    if isfield(S,'TMT_linear') && ~isempty(S.TMT_linear) && ~S.TMT_linear
        if S.tmt && isfield(S,'k1TMT') && ~isempty(S.k1TMT) && ~S.tmt_locked
            savenameparts{end+1} = ['k' num2str(S.k1TMT)];
            casfuncfolparts{end+1} = ['k' num2str(S.k1TMT)];
        end
        if S.tmt && isfield(S,'k2TMT') && ~isempty(S.k2TMT) && ~S.tmt_locked
            savenameparts{end+1} = ['kc' num2str(S.k2TMT)];
            casfuncfolparts{end+1} = ['kc' num2str(S.k2TMT)];
        end
        if S.tmt && isfield(S,'t1TMT') && ~isempty(S.t1TMT) && ~S.tmt_locked
            savenameparts{end+1} = ['t' num2str(S.t1TMT*10)];
            casfuncfolparts{end+1} = ['t' num2str(S.t1TMT*10)];
        end
    else
        if S.tmt && isfield(S,'kTMT') && ~isempty(S.kTMT) && ~S.tmt_locked
            savenameparts{end+1} = ['k' num2str(S.kTMT)];
            casfuncfolparts{end+1} = ['k' num2str(S.kTMT)];
        end
    end
    if isfield(S,'Windlass') && ~isempty(S.Windlass) 
        if S.Windlass
            savenameparts{end+1} = ['WL'];
            casfuncfolparts{end+1} = ['WL'];
            if isfield(S,'cWL') && ~isempty(S.cWL)
                savenameparts{end} = [savenameparts{end} num2str(S.cWL*1000)];
                casfuncfolparts{end} = [casfuncfolparts{end} num2str(S.cWL*1000)];
            end
        else
            not{end+1} = 'not_WL';
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
if isfield(S,'ExoBool') && ~isempty(S.ExoBool)
    if S.ExoBool == 1
        if isfield(S,'ExoScale') && ~isempty(S.ExoScale)
            if S.ExoScale == 0
                savenameparts{end+1} = ['pas'];
            else
                savenameparts{end+1} = ['act'];
                if isfield(S,'ExoImplementation') && ~isempty(S.ExoImplementation)
                    if strcmp(S.ExoImplementation,'TorqueTibiaCalcn')
                        savenameparts{end+1} = ['TTC'];
                    elseif strcmp(S.ExoImplementation,'TorqueTibiaMetatarsi')
                        savenameparts{end+1} = ['TTM'];
                    end
                end
                if isfield(S,'ExoController') && ~isempty(S.ExoController) 
                    if strcmp(S.ExoController,'Ideal Assistance')
                        savenameparts{end+1} = ['IA'];
                        if isfield(S,'T_max_ankle_exo') && ~isempty(S.T_max_ankle_exo)
                            savenameparts{end} = [savenameparts{end} 'Tx' num2str(S.T_max_ankle_exo)];
                        end
                        if isfield(S,'T_min_ankle_exo') && ~isempty(S.T_max_ankle_exo)
                            savenameparts{end} = [savenameparts{end} 'Tn' num2str(S.T_min_ankle_exo)];
                        end
                        if isfield(S,'P_max_ankle_exo') && ~isempty(S.P_max_ankle_exo)
                            savenameparts{end} = [savenameparts{end} 'P' num2str(S.P_max_ankle_exo)];
                        end
                    else
                        not{end+1} = 'not_IA';
                    end
                end
            end
        end
    else
        not{end+1} = 'not_pas';
        not{end+1} = 'not_act';
    end
end

savename = savenameparts{1};
for i=2:numel(savenameparts)
    savename = [savename '_' savenameparts{i}];
end
criteria = {};
for i=1:numel(savenameparts)
    if ~(~isfield(S,'MuscModelAsmp') && i==nr)
        criteria{end+1} = savenameparts{i};
    end
end
casfuncfol = casfuncfolparts{1};
for i=2:numel(casfuncfolparts)
    casfuncfol = [casfuncfol '_' casfuncfolparts{i}];
end

for i=1:numel(not)
    criteria{end+1} = not{i};
end

% disp(savename);
% disp(casfuncfol);


if nargout == 3
    varargout{1} = criteria;
end

end