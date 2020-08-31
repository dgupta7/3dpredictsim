function [massM,tensions] = GetMuscleMass(muscleNames,params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% get the muscle parameters
FMo = params(1,:);
lMo = params(2,:);

% get the specific tension for each muscle
tension = getSpecificTensions(muscleNames(1:end-3));
tensions = [tension;tension]';

% compute the mass of the muscle
volM = FMo.*lMo;
massM = volM.*(1059.7)./(tensions*1e6);

end

