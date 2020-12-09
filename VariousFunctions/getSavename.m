function [savename, casfuncfol] = getSavename(S)
% build savename and casadifunction foldername
if strcmp(S.subject,'s1_Poggensee')
    savename = 'Pog_s1';
    casfuncfol = 'casadi_s1Pog';
end
if S.tmt
    savename = [savename '_tmt'];
    casfuncfol = [casfuncfol '_tmt'];
    if S.tmt_locked
        savename = [savename 'L'];
    end
end
if isfield(S,'MuscModelAsmp') && ~isempty(S.MuscModelAsmp) && S.MuscModelAsmp==0
    savename = [savename '_bCst'];
    casfuncfol = [casfuncfol '_MuscModel_bCst'];
else
    savename = [savename '_aCst'];
    casfuncfol = [casfuncfol '_MuscModel_alphaCst'];
end
if S.tmt && isfield(S,'dTMT') && ~isempty(S.dTMT) && ~S.tmt_locked
    savename = [savename '_d0' num2str(S.dTMT*10)];
    casfuncfol = [casfuncfol '_d0' num2str(S.dTMT*10)];
end
if S.tmt && isfield(S,'kTMT') && ~isempty(S.kTMT) && ~S.tmt_locked
    savename = [savename '_k' num2str(S.kTMT)];
    casfuncfol = [casfuncfol '_k' num2str(S.kTMT)];
end
if S.IGsel == 1
    savename = [savename '_ig1'];
else
    savename = [savename '_ig2' num2str(S.IGmodeID )];
end
if S.ExoBool == 1
    if S.ExoScale == 0
        savename = [savename '_pas'];
    else
        savename = [savename '_act'];
    end    
end
disp(savename);
disp(casfuncfol);
end