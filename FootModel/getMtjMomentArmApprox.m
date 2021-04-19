function [MA_mtj] = getMtjMomentArmApprox(q_mtj)
% Returns approximated momentarms for the muscles spanning the mt joint,
% around the mt joint.
%
% Author: Lars D'Hondt (April 2021)


%% flex_dig
% run \FootModel\solveFootmodelParameters.m to get these values
flex_dig2mtj = 0.042381;
mtj2flex_dig = 0.0092642;
phi0 = 0.032319;

phi = phi0 + q_mtj;

l = sqrt(flex_dig2mtj^2+mtj2flex_dig^2-2*flex_dig2mtj*mtj2flex_dig.*cos(phi));
MA = flex_dig2mtj*mtj2flex_dig./l.*sin(phi);

MA_1 = -MA;

%% flex_hal
% run \FootModel\solveFootmodelParameters.m to get these values
flex_hal2mtj = 0.037305;
mtj2flex_hal = 0.031181;
phi0 = 2.5322;

phi = phi0 + q_mtj;

l = sqrt(flex_hal2mtj^2+mtj2flex_hal^2-2*flex_hal2mtj*mtj2flex_hal.*cos(phi));
MA = flex_hal2mtj*mtj2flex_hal./l.*sin(phi);

MA_2 = -MA;

%% ext_dig
ext_dig2mtj = 0.027608;

MA_3 = ext_dig2mtj;

%% ext_hal
ext_hal2mtj = 0.037684;

MA_4 = ext_hal2mtj;



%%
MA_mtj = [MA_1; MA_2; MA_3; MA_4];


end