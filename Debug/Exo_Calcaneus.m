

%% Exoskeleton actuation at tibia and calcaneus
%----------------------------------------------


% Here we want to implement the actuation of the exoskeleton at the level
% of the calcaneus and the tibia (instead of a torque at the ankle joint)

% We can do this by applying a torque of equal magnitude, but in opposite
% direction at the level of the tibia and calacneus. However, this is
% harder to implement in c++, especially since the axis of the subtalar
% joint is not aligned with the ankle axis

% Notes on Anatomy:
% Tibia - (Ankle) - Talus - (subtalar) - calcaneus

% I believe that in this case, we can do the implementation in matlab. We
% add to the outcome of ID an ankle torque to the ankle axis. We also add a
% torque to the subtalar axis, which is equal to the ankle torque projected
% on the subtalar joint.


AxisAnkle = [ 0.10501355, 0.17402245, 0.97912632];  % orientation ankle axis in frame of tibia
AxisSubt  = [ -0.78717961, -0.60474746, -0.12094949]; % orientation subtalar axis in frame of talus
% AxisSubt  = [0 0 1]; % orientation subtalar axis in frame of talus
Texo = [0 0 -10];
Tsubt = Texo*AxisSubt';
disp(Tsubt);





