function [M_mtj, M_mtpj, varargout] = getPassiveMtjMomentWindlass_v3(q_mt,qdot_mt,q_mtp,f_PF_stiffness,varargin)
% This function returns total the passive moment around the tarsometatarsal
% joint, and the extra passive moment around the mtp joint caused be the
% plantar fascia.
% Input arguments:
%   q_mt: angle of midtarsal joint (rad)
%   q_mtp: angle of mtp joint(rad)
%   qdot_mt: angular velocity of midtarsal joint (rad/s)
%   f_PF_stiffness: function that takes the plantar fascia length, and
%                   returns the force.
%   varargin: possibility to pass the settings struct S
% 
% Author: Lars D'Hondt (March 2021)

%% Geometry
% run \FootModel\solveFootmodelParameters.m to get these values
% foot arch
calcn2mtj = 0.08207;
mtj2mtpj = 0.089638;
beta0 = 2.4935;
% windlass mechanism
calcnPF2mtj = 0.06695;
mtj2mttPF = 0.091714;
phi0 = 2.1274;


%% Default parameters
R_mtth = 7.5e-3; % radius of the metatarsal head
l_toe = 0; % distance from metatarsal head to PF attachment point at toe
nl = 1; % use nonlinear ligament stiffness
kMT_li = 90; % linear stiffness

%% Get specific parameters
if isempty(varargin)
    dMT = 0; 
else
    dMT = S.dTMT;
    nl = S.MT_li_nonl;
    kMT_li = S.kMT_li;
end


%% Geometry Windlass mechanism
% based on midtarsal joint angle
phi = phi0 + q_mt; % top angle of WL triangle
l_PF_fa = sqrt(calcnPF2mtj^2 + mtj2mttPF^2 - 2*calcnPF2mtj*mtj2mttPF*cos(phi)); % length of PF spanning arch
MA_PF = calcnPF2mtj*mtj2mttPF/l_PF_fa*sin(phi); % moment arm of PF to mtj

% Plantar fascia length
l_PF = l_PF_fa + R_mtth*q_mtp + l_toe;

%% Torques
% plantar fascia
F_PF = f_PF_stiffness(l_PF);
if nargout > 1
    % convert SX to double when called from post-processing
    F_PF = full(F_PF);
end

M_PF = -F_PF*MA_PF;

% other elastic structures
if nl
    k1 = 90;
    dl_0 = 3*pi/180;
    dq_0 = 2*pi/180;
    M_li = -k1*( (q_mt-dq_0) - dl_0*tanh((q_mt-dq_0)/dl_0/1.2)) +... % toe-in and linear part
            (-exp(20*(q_mt-20*pi/180)) + exp(-25*(q_mt+10*pi/180)))/2; % stiffening at the end
    
else
    M_li = -kMT_li*q_mt;
end

% viscous damping
M_d = -dMT*qdot_mt;

% Total passive moment
M_mtj = M_PF + M_li + M_d;

% Reaction torque on toes
M_mtpj = -R_mtth*F_PF;

%% Return some more outputs if the function is called from postprocessing
if nargout > 2
    % Geometry foot arch
    % based on midtarsal joint angle
    beta = beta0 + q_mt; % top angle of foot arch triangle
    l_fa = sqrt(calcn2mtj^2 + mtj2mtpj^2 - 2*calcn2mtj*mtj2mtpj*cos(beta)); % foot arch length
    h_fa = calcn2mtj*mtj2mtpj/l_fa*sin(beta); % foot arch height
    
    varargout = {M_PF,M_li,M_d,F_PF,l_PF,MA_PF,l_PF_fa,h_fa,l_fa};
end






end