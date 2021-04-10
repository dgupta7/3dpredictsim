function [MA_mtp] = getMtpMomentArmApprox(q_mtpj)
% Returns approximated momentarms for the muscles spanning the mtp joint,
% around the mtp joint.
%
% Author: Lars D'Hondt (April 2021)


%% flex_dig
MA_1 = -0.007;

%% flex_hal
MA_2 = -0.006;

%% ext_dig
% run \FootModel\solveFootmodelParameters.m to get these values
ext_dig2mtpj = 0.0091016;
mtpj2ext_dig = 0.028891;
phi0 = 0.59134;

phi = phi0 + q_mtpj;
l = sqrt(ext_dig2mtpj^2 + mtpj2ext_dig^2 - 2*ext_dig2mtpj*mtpj2ext_dig*cos(phi));
MA_3 = ext_dig2mtpj*mtpj2ext_dig/l*sin(phi);


%% ext_hal
% run \FootModel\solveFootmodelParameters.m to get these values
ext_hal2mtpj = 0.017156;
mtpj2ext_hal = 0.033969;
phi0 = 0.59189;

phi = phi0 + q_mtpj;
l = sqrt(ext_hal2mtpj^2 + mtpj2ext_hal^2 - 2*ext_hal2mtpj*mtpj2ext_hal*cos(phi));
MA_4 = ext_hal2mtpj*mtpj2ext_hal/l*sin(phi);


MA_mtp = [MA_1; MA_2; MA_3; MA_4];


end