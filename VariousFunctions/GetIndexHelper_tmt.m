function [IndexLeft,IndexRight,QsInvA,QsInvB,QdotsInvA...
    ,QdotsInvB,orderQsOpp] = GetIndexHelper_tmt(S,jointi)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

if strcmp(S.ModelName,'Rajagopal')
    % indexes to select kinematics left and right leg
    IndexLeft = [jointi.hip_flex.l jointi.hip_add.l jointi.hip_rot.l, ...
        jointi.knee.l jointi.ankle.l jointi.subt.l jointi.mtp.l];
    IndexRight = [jointi.hip_flex.r jointi.hip_add.r jointi.hip_rot.r, ...
        jointi.knee.r jointi.ankle.r jointi.subt.r jointi.mtp.r];
elseif strcmp(S.ModelName,'Gait92')
    if isfield(S,'FootMuscles') && S.FootMuscles
        % indexes to select kinematics left and right leg
        IndexLeft = [jointi.hip_flex.l jointi.hip_add.l jointi.hip_rot.l, ...
            jointi.knee.l jointi.ankle.l jointi.subt.l jointi.mtj.l jointi.mtp.l,...
            jointi.trunk.ext, jointi.trunk.ben, jointi.trunk.rot];
        IndexRight = [jointi.hip_flex.r jointi.hip_add.r jointi.hip_rot.r, ...
            jointi.knee.r jointi.ankle.r jointi.subt.r jointi.mtj.r jointi.mtp.r,...
            jointi.trunk.ext, jointi.trunk.ben, jointi.trunk.rot];

    else
        % indexes to select kinematics left and right leg
        IndexLeft = [jointi.hip_flex.l jointi.hip_add.l jointi.hip_rot.l, ...
            jointi.knee.l jointi.ankle.l jointi.subt.l jointi.mtp.l,...
            jointi.trunk.ext, jointi.trunk.ben, jointi.trunk.rot];
        IndexRight = [jointi.hip_flex.r jointi.hip_add.r jointi.hip_rot.r, ...
            jointi.knee.r jointi.ankle.r jointi.subt.r jointi.mtp.r,...
            jointi.trunk.ext, jointi.trunk.ben, jointi.trunk.rot];
    end
end

% indexes for symmetry steps
QsInvA = [jointi.pelvis.tilt,...
    jointi.pelvis.ty,...
    jointi.hip_flex.l:jointi.trunk.ext,...
    jointi.sh_flex.l:jointi.elb.r]';
QsInvB = [jointi.pelvis.tilt,...
    jointi.pelvis.ty,...
    jointi.hip_flex.r:jointi.hip_rot.r,...
    jointi.hip_flex.l:jointi.hip_rot.l,...
    jointi.knee.r,...
    jointi.knee.l,...
    jointi.ankle.r,...
    jointi.ankle.l,...
    jointi.subt.r,...
    jointi.subt.l,...
    jointi.tmt.r,...
    jointi.tmt.l,...
    jointi.mtp.r,...
    jointi.mtp.l,...
    jointi.trunk.ext,...
    jointi.sh_flex.r:jointi.sh_rot.r,...
    jointi.sh_flex.l:jointi.sh_rot.l,...
    jointi.elb.r,...
    jointi.elb.l]';

QdotsInvA = [jointi.pelvis.tilt,...
    jointi.pelvis.tx,jointi.pelvis.ty,...
    jointi.hip_flex.l:jointi.trunk.ext,...
    jointi.sh_flex.l:jointi.elb.r]';
QdotsInvB = [jointi.pelvis.tilt,...
    jointi.pelvis.tx,jointi.pelvis.ty,...
    jointi.hip_flex.r:jointi.hip_rot.r,...
    jointi.hip_flex.l:jointi.hip_rot.l,...
    jointi.knee.r,...
    jointi.knee.l,...
    jointi.ankle.r,...
    jointi.ankle.l,...
    jointi.subt.r,...
    jointi.subt.l,...
    jointi.tmt.r,...
    jointi.tmt.l,...
    jointi.mtp.r,...
    jointi.mtp.l,...
    jointi.trunk.ext,...
    jointi.sh_flex.r:jointi.sh_rot.r,...
    jointi.sh_flex.l:jointi.sh_rot.l,...
    jointi.elb.r,...
    jointi.elb.l]';

orderQsOpp = [jointi.pelvis.list:jointi.pelvis.list,...
    jointi.pelvis.rot:jointi.pelvis.rot,...
    jointi.pelvis.tz:jointi.pelvis.tz,...
    jointi.trunk.ben:jointi.trunk.ben,...
    jointi.trunk.rot:jointi.trunk.rot];



end

