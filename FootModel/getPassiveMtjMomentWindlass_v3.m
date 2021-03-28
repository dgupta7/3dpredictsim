function [M_mtj, M_mtpj, varargout] = getPassiveMtjMomentWindlass_v3(q_mt,qdot_mt,q_mtp,varargin)
% This function returns the passive moment around the tarsometatarsal joint
% q_mt: angle of midtarsal joint (rad)
% q_mtp: angle of mtp joint(rad)
% qdot_mt: angular velocity of midtarsal joint (rad/s)
% kMT_li: angular stiffness of mt joint, from other ligaments (Nm/rad)
% dMT: angular viscous friction of tmt joint (Nms/rad)
% 
% Author: Lars D'Hondt (March 2021)

%% Geometry
% run \FootModel\solveFootmodelParameters.m to get these values
% windlass mechanism
calcn2mtj = 0.08207;
mtj2mtpj = 0.089638;
beta0 = 2.4935;
% foot arch
calcnPF2mtj = 0.06695;
mtj2mttPF = 0.091714;
phi0 = 2.1274;


%% Default parameters

R_mtth = 7.5e-3; % radius of the metatarsal head
l_toe = 0; % distance from metatarsal head to PF attachment point at toe

% The reference position of the mtj (used for calculating the torque) is
% given by: a + b * q_mtp
q_mt_0_a = 5*pi/180;
q_mt_0_b = -0;


%% Get specific parameters
if isempty(varargin)
    f_PF_stiffness = @(le) 5.9706e+05*(le-0.17); % linear with default values
else
    f_PF_stiffness = varargin{1};
end

kMT_li = 0;
dMT = 0;


if length(varargin)==5
    % when passed separately
    kMT_li = varargin{2};
    dMT = varargin{3};
    
elseif length(varargin)==2
    % from struct with settings
    S = varargin{2};
    kMT_li = S.kTMT_li;
    dMT = S.dTMT;
    
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
q0 = q_mt_0_a + q_mt_0_b*q_mtp; % neutral position
M_li = -(q_mt - q0)*kMT_li; % elastic moment

% % joint stiffens close to limit positions
% k_pass.mtj = [-1 20 1 -20]';
% theta.pass.mtj = [-0.20 0.20]';
% tau_pass = k_pass.mtj(1,1)*exp(k_pass.mtj(2,1)*(q_mt-theta.pass.mtj(2,1))) + ...
%     k_pass.mtj(3,1)*exp(k_pass.mtj(4,1)*(q_mt-theta.pass.mtj(1,1)));
% M_li = M_li + tau_pass;

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