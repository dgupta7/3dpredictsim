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

if isfield(S,'tanh_b') && ~isempty(S.tanh_b)
    savenameparts{end+1} = ['tanh' num2str(S.tanh_b)];
end

if isfield(S,'mtj') && ~isempty(S.mtj) && S.mtj
    if isfield(S,'PF_stiffness') && ~isempty(S.PF_stiffness)
        savenameparts{end+1} = ['PF_' S.PF_stiffness];
        casfuncfolparts{end+1} = ['PF_' S.PF_stiffness];
    end
    if isfield(S,'PF_slack_length') && ~isempty(S.PF_slack_length)
        savenameparts{end+1} = ['ls' num2str(S.PF_slack_length*1000)];
        casfuncfolparts{end+1} = ['ls' num2str(S.PF_slack_length*1000)];
    end
    if isfield(S,'MT_li_nonl') && ~isempty(S.MT_li_nonl) && S.MT_li_nonl
        if isfield(S,'mtj_stiffness') && ~isempty(S.mtj_stiffness)
            if strcmp(S.mtj_stiffness,'signed_lin')
                savenameparts{end+1} = ['MT_nl_k' num2str(S.kMT_li) '_' num2str(S.kMT_li2)];
                casfuncfolparts{end+1} = ['MT_nl_k' num2str(S.kMT_li) '_' num2str(S.kMT_li2)];
            else
                savenameparts{end+1} = ['MT_nl_' S.mtj_stiffness];
                casfuncfolparts{end+1} = ['MT_nl_' S.mtj_stiffness];
            end
        else
            savenameparts{end+1} = ['MT_nl'];
            casfuncfolparts{end+1} = ['MT_nl'];
        end
    elseif isfield(S,'kMT_li') && ~isempty(S.kMT_li)
        savenameparts{end+1} = ['MT_k' num2str(S.kMT_li)];
        casfuncfolparts{end+1} = ['MT_k' num2str(S.kMT_li)];
        
    end
    if isfield(S,'MT_li_nonl') && ~isempty(S.MT_li_nonl) && ~S.MT_li_nonl
        not{end+1} = 'not_MT_nl';
    end
    if isfield(S,'dMT') && ~isempty(S.dMT) && S.dMT~=0
        savenameparts{end+1} = ['d0' num2str(S.dMT*10)];
        casfuncfolparts{end+1} = ['d0' num2str(S.dMT*10)];
    end
    if isfield(S,'WL_T_mtp') && ~isempty(S.WL_T_mtp)
        savenameparts{end+1} = 'MTP';
        casfuncfolparts{end+1} = 'MTP';
        if S.WL_T_mtp
            if isfield(S,'FootMuscles') && ~isempty(S.FootMuscles) && S.FootMuscles
                savenameparts{end} = [savenameparts{end} '_Mf'];
                casfuncfolparts{end} = [casfuncfolparts{end} '_Mf'];
            elseif isfield(S,'Mu_mtp') && ~isempty(S.Mu_mtp) && S.Mu_mtp
                savenameparts{end} = [savenameparts{end} '_Mu'];
                casfuncfolparts{end} = [casfuncfolparts{end} '_Mu'];
            else
                savenameparts{end} = [savenameparts{end} '_T'];
                casfuncfolparts{end} = [casfuncfolparts{end} '_T'];
            end
        else
            savenameparts{end} = [savenameparts{end} '_k'];
            casfuncfolparts{end} = [casfuncfolparts{end} '_k'];
        end
        if isfield(S,'kMTP') && ~isempty(S.kMTP)
            savenameparts{end} = [savenameparts{end} num2str(S.kMTP)];
            casfuncfolparts{end} = [casfuncfolparts{end} num2str(S.kMTP)];
        end
        if isfield(S,'dMTP') && ~isempty(S.dMTP)
        savenameparts{end+1} = ['d0' num2str(S.dMTP*10)];
        casfuncfolparts{end+1} = ['d0' num2str(S.dMTP*10)];
    end
    end
    if isfield(S,'stiffen_arch') && ~isempty(S.stiffen_arch) && S.stiffen_arch~=0
        savenameparts{end+1} = ['K' num2str(S.stiffen_arch)];
        casfuncfolparts{end+1} = ['K' num2str(S.stiffen_arch)];
    end
    if isfield(S,'contactStiff') && ~isempty(S.contactStiff)
        if S.contactStiff == 10
            savenameparts{end+1} = 'spx10';
        elseif S.contactStiff == 2
            savenameparts{end+1} = 'spx2';
        end
    end
    if isfield(S,'WLpoly') && ~isempty(S.WLpoly) && S.WLpoly==0
        savenameparts{end+1} = 'np';
        casfuncfolparts{end+1} = 'np';
    end
    if isfield(S,'PIM') && ~isempty(S.PIM) && S.PIM==1
        savenameparts{end+1} = 'PIM';
        casfuncfolparts{end+1} = 'PIM';
        if isfield(S.W,'PIM') && ~isempty(S.W.PIM)
            savenameparts{end+1} = ['w' num2str(S.W.PIM,2)];
        end
        if isfield(S.W,'P_PIM') && ~isempty(S.W.P_PIM)
            savenameparts{end+1} = ['w' num2str(S.W.P_PIM,2)];
        end
    end
else
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
        if isfield(S,'Windlass') && ~isempty(S.Windlass) && S.tmt
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
                    elseif strcmp(S.ExoImplementation,'TorqueTibiaCalcnMetatarsi')
                        savenameparts{end+1} = ['TTCM'];
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