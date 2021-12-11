% This function retrieves the length and moment arms of a ligament from an
% OpenSim model. The ligament is hard coded to be the right side plantar
% fascia.
%
% Author: Lars D'Hondt
% Date: 7/Dec/2021
%
function [] = getPlantarFasciaGeometry(MainPath,Modelpath,PolyFolder)

% Output folder for this subject
OutFolder = PolyFolder;
SubjFolder = fullfile(MainPath,'Polynomials',OutFolder);

%% Initialise model
import org.opensim.modeling.*;
model = Model(Modelpath);
s = model.initSystem;

%% Get plantar fascia length in neutral position
% Set all state variables to 0
state_vars = model.getStateVariableValues(s);
state_vars.setToZero();
model.setStateVariableValues(s,state_vars);
model.realizePosition(s);
% Get plantar fascia length
PF = Ligament.safeDownCast(model.getForceSet().get('PlantarFascia_r'));
l_0 = PF.getLength(s);


%% Find coordinates that interact with the plantar fascia
idx_PF = []; % coordinate indices in the OpenSim API
dof_crossing = {}; % coordinate names
% Loop through coordinates
for i=0:model.getCoordinateSet().getSize()-1
    % Set the positional state of this coordinate to 5°
    state_vars.setToZero();
    state_vars.set(2*i,5*pi/180);
    model.setStateVariableValues(s,state_vars);
    model.realizePosition(s);
    % If the plantar fascia length changed, it spans this coordinate
    l_check = PF.getLength(s);
    if abs(l_check-l_0) > 1e-6
        idx_PF(end+1) = i;
        dof_crossing{end+1} = char(model.getCoordinateSet().get(i).getName());
    end
end

%% Create dummy motion
% Boundaries on dofs, format: name, min, max (°)
bounds = {'mtp_angle_r',-30,50;
        'mtj_angle_r',-20,20;
        'tmt_angle_r',-10,10};
% Random table
n = 500;
p = length(idx_PF);
X = lhsdesign(n,p);
% Apply bounds
X_scale = zeros(1,p);
X_min = zeros(1,p);
for i=1:p
    idx_bounds = find(strcmp(bounds(:,1),dof_crossing{i}));
    X_scale(i) = abs(bounds{idx_bounds,2}-bounds{idx_bounds,3});
    X_min(i) = bounds{idx_bounds,2};
end
% Coordinate angles (in rad)
Angles = X.*(ones(n,1)*X_scale)+(ones(n,1)*X_min);
Angles = Angles*pi/180;

%% Evaluate plantar fascia lenght and moment arms for dummy motion
% Set state vector to 0
state_vars = model.getStateVariableValues(s);
state_vars.setToZero();
model.setStateVariableValues(s,state_vars);
model.realizePosition(s);
% Initialise matrices for results
l_PF = zeros(n,1);
d_PF = zeros(n,1,p);
% Loop through dummy states
for i=1:n
    % Set each coordinate value
    for j=1:p
        state_vars.set(idx_PF(j)*2,Angles(i,j));
    end
    model.setStateVariableValues(s,state_vars);
    model.realizePosition(s);
    % Get plantar fascia length
    l_PF(i,:) = PF.getLength(s);
    % Get moment arm for each joint
    for j=1:p
        d_PF(i,1,j) = PF.computeMomentArm(s,model.getCoordinateSet().get(idx_PF(j)));
    end

end

%% Combine results in struct
% Field names are the same as the data resulting from muscle analysis, so
% they can be processed by the same function.
LigamentData.dof_names = dof_crossing;
LigamentData.muscle_names = {'Plantar_Fascia_r'};
LigamentData.lMT = l_PF;
LigamentData.dM = d_PF;
LigamentData.q = Angles;

%% Call PolynomialFit
    [ligament_spanning_joint_INFO,LigamentInfo] = PolynomialFit_mtp(LigamentData);
    save(fullfile(SubjFolder,'LigamentData.mat'),'LigamentData')
    save(fullfile(SubjFolder,'ligament_spanning_joint_INFO.mat'),'ligament_spanning_joint_INFO')
    save(fullfile(SubjFolder,'LigamentInfo.mat'),'LigamentInfo');








