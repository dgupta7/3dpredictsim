function [M, varargout] = getPassiveMtjMomentWindlass_v2(q_mt,qdot_mt,q_mtp,varargin)
% This function returns the passive moment around the tarsometatarsal joint
% q_mt: angle of midtarsal joint (rad)
% q_mtp: angle of mtp joint(rad)
% qdot_mt: angular velocity of midtarsal joint (rad/s)
% kMT_li: angular stiffness of mt joint, from other ligaments (Nm/rad)
% kMT_PA: angular stiffness of mt joint, from plantar fascia (Nm/rad)
%           Formulated like this for consistency with earlier models
% dMT: angular viscous friction of tmt joint (Nms/rad)
% Subject: test subject
% cWL: relation between mtp angle and foot arch length
%           (see doi.org/10.1098/rsif.2018.0270)

%% Get default parameters
a = 0.08207;
b = 0.089638;
phi0 = 2.493499;
H0 = 0.027280;

% a = 0.0857;
% b = 0.0932;
% phi0 = 2.2799;
% H0 = 0.0373;

R_mtth = 0.01; % radius of the metatarsal head
l_toe = 0.01; % distance from metatarsal head to PF attachment point at toe

% The reference position of the mtj (used for calculating the torque) is
% given by: a + b * q_mtp
q_mt_0_a = 0*pi/180;
q_mt_0_b = -0;


%% Get specific parameters
if isempty(varargin)
    f_PF_stiffness = @(le) 5.9706e+05*(le-0.17); % linear with default values
else
    f_PF_stiffness = varargin{1};
end

kMT_li = 0;
dMT = 0;
cWL = 0.2;

if length(varargin)==5
    % when passed separately
    kMT_li = varargin{2};
    dMT = varargin{3};
    subject = varargin{4};
    cWL = varargin{5};
    
elseif length(varargin)==2
    % from struct with settings
    S = varargin{2};
    kMT_li = S.kTMT_li;
    dMT = S.dTMT;
    subject = S.subject;
    cWL = S.cWL;
end

% Get subject-specific constants, derived from foot geometry (projected on
% sagittal plane)
% (see \FootModel\WindlassParameters.m)
% if strcmp(subject,'s1_Poggensee')
%     H0 = 0.0724; % unloaded foot arch height (m)
%     a = 0.1126; % length from calcn origin to mtj (m)
%     b = 0.1047; % length from mtj to mtpj (m)
%     phi0 = 1.6798; % angle between a and b (rad)
% elseif strcmp(subject,'subject1')
%     H0 = 0.0724; % unloaded foot arch height (m)
%     a = 0.1126;
%     b = 0.1047;
%     phi0 = 1.6798;
% end

%% Geometry Windlass mechanism
% zero load
L0_fa = sqrt(a^2 + b^2 - 2*a*b*cos(phi0)); % foot arch length for mtp angle 0
l_0_fa = (1-cWL*(q_mtp*180/pi)/20)*L0_fa; % foot arch length
h_0_fa = a*b./l_0_fa*sin(phi0);
q_mt_0 = acos( (a^2 + b^2 - l_0_fa.^2)/(2*a*b) ) - phi0;

% based on midtarsal joint angle
phi = phi0 + q_mt;
l_fa = sqrt(a^2 + b^2 - 2*a*b*cos(phi));
h = a*b./l_fa.*sin(phi);

% Plantar fascia length
l = l_fa + R_mtth*q_mtp + l_toe;

%% Torques
% plantar fascia
if nargout == 1
    F_PF = f_PF_stiffness(l); % when building passive torques function
else
    F_PF = full(f_PF_stiffness(l)); % when post-processing
end
M_PF = F_PF*h;

% other elastic structures
q0 = q_mt_0_a + q_mt_0_b*q_mtp;
M_li = (q_mt - q0)*kMT_li;

% viscous damping
M_d = dMT*qdot_mt;

%% Total passive moment
M = -(M_PF + M_li + M_d);

% Return some more outputs if the function is called from postprocessing
if nargout == 2
    varargout = {q_mt_0};
elseif nargout > 2
    varargout = {M_PF,F_PF,M_li,M_d,l,l_fa,l_0_fa,L0_fa,h,h_0_fa,H0,q_mt_0};
end

end