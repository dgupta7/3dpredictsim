function [] = compare_Lundgren_2008(ResultsFile)

load(ResultsFile,'R');


Modelpath = 'C:\Users\u0150099\Documents\master_thesis\3dpredictsim\OpenSimModel\subject1\Fal_s1_mtj_sc.osim';

istance = 1:1:ceil(R.Event.Stance);
x = linspace(1,100,length(istance));

q_ankle = R.Qs(istance,strcmp(R.colheaders.joints,'ankle_angle_r'))*pi/180;
q_subt = R.Qs(istance,strcmp(R.colheaders.joints,'subtalar_angle_r'))*pi/180;
q_mtj = R.Qs(istance,strcmp(R.colheaders.joints,'mtj_angle_r'))*pi/180;

%% Initialise model
import org.opensim.modeling.*;
model = Model(Modelpath);
s = model.initSystem;

coord_names = {model.getCoordinateSet().getSize()};
for i=1:model.getCoordinateSet().getSize()
    coord_names{i} = char(model.getCoordinateSet().get(i-1));
end

iankle = find(strcmp(coord_names,'ankle_angle_r'))-1;
isubt = find(strcmp(coord_names,'subtalar_angle_r'))-1;
imtj = find(strcmp(coord_names,'mtj_angle_r'))-1;

% Get bodies
tibia = model.getBodySet().get('tibia_r');
talus = model.getBodySet().get('talus_r');
calcn = model.getBodySet().get('calcn_r');
midfoot = model.getBodySet().get('midfoot_r');
forefoot = model.getBodySet().get('forefoot_r');

% Set state vector to 0
state_vars = model.getStateVariableValues(s);
state_vars.setToZero();
model.setStateVariableValues(s,state_vars);
model.realizePosition(s);

%% loop through motion
% Angles given as: [adduction, evertion, dorsiflexion]

talus_rt_tibia = zeros(length(q_ankle),3);
calcn_rt_tibia = zeros(length(q_ankle),3);
calcn_rt_talus = zeros(length(q_ankle),3);
mf_rt_talus = zeros(length(q_ankle),3);
mf_rt_calcn = zeros(length(q_ankle),3);
ff_rt_talus = zeros(length(q_ankle),3);

for i=1:length(q_ankle)
    % Set each coordinate value
    state_vars.set(iankle*2,q_ankle(i));
    state_vars.set(isubt*2,q_subt(i));
    state_vars.set(imtj*2,q_mtj(i));
    model.setStateVariableValues(s,state_vars);
    model.realizePosition(s);
    

    talus_rt_tibia(i,:) = talus.findTransformBetween(s,tibia).R().convertRotationToBodyFixedXYZ().getAsMat'*180/pi;
    calcn_rt_tibia(i,:) = calcn.findTransformBetween(s,tibia).R().convertRotationToBodyFixedXYZ().getAsMat'*180/pi;
    calcn_rt_talus(i,:) = calcn.findTransformBetween(s,talus).R().convertRotationToBodyFixedXYZ().getAsMat'*180/pi;
    mf_rt_talus(i,:) = midfoot.findTransformBetween(s,talus).R().convertRotationToBodyFixedXYZ().getAsMat'*180/pi;
    mf_rt_calcn(i,:) = midfoot.findTransformBetween(s,calcn).R().convertRotationToBodyFixedXYZ().getAsMat'*180/pi;
    ff_rt_talus(i,:) = forefoot.findTransformBetween(s,talus).R().convertRotationToBodyFixedXYZ().getAsMat'*180/pi;

end



%%

figure

subplot(6,3,1)
plot(x,-talus_rt_tibia(:,3))
hold on
ylabel('-DF/+PF (°)')
title('talus relative to tibia')

subplot(6,3,4)
plot(x,talus_rt_tibia(:,2))
hold on
ylabel('-INV/+EV (°)')

subplot(6,3,7)
plot(x,talus_rt_tibia(:,1))
hold on
ylabel('-ABD/+ADD (°)')

%
subplot(6,3,2)
plot(x,-calcn_rt_tibia(:,3))
hold on
ylabel('-DF/+PF (°)')
title('calcaneus relative to tibia')

subplot(6,3,5)
plot(x,calcn_rt_tibia(:,2))
hold on
ylabel('-INV/+EV (°)')

subplot(6,3,8)
plot(x,calcn_rt_tibia(:,1))
hold on
ylabel('-ABD/+ADD (°)')


%
subplot(6,3,3)
plot(x,-calcn_rt_talus(:,3))
hold on
ylabel('-DF/+PF (°)')
title('calcaneus relative to talus')

subplot(6,3,6)
plot(x,calcn_rt_talus(:,2))
hold on
ylabel('-INV/+EV (°)')

subplot(6,3,9)
plot(x,calcn_rt_talus(:,1))
hold on
ylabel('-ABD/+ADD (°)')















