%% Convert units

sInch = 0.0254;
slbs  = 0.45359237;

% convert rotational inertia -- Exo Foot
IFoot = [7.32; 17.02; 23.25] ;
iFoot_si = IFoot.*(slbs.*sInch^2)

% convert rotational inertia -- Exo Tibia
ITibia = [24.91; 22.69; 9.32] ;
iTibia_si = ITibia.*(slbs.*sInch^2)

% axes information - Poggensee
% x= forward
% y = to the left
% z = up

% axes information opensim
% x = forward
% y = upward
% z = right

% (difference between both is 90° rotation around x-axis

iFoot_osim = abs(rotx(90)*iFoot_si)
iTibia_osim = abs(rotx(90)*iTibia_si)



