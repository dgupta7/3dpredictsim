function [M, varargout] = getPassiveMtjMomentWindlass(q_mt,qdot_mt,q_mtp,kTMT_li,kMT_PF,dMT,subject,cWL)
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


% Get subject-specific constants, derived from foot geometry (projected on
% sagittal plane)
% (see \FootModel\WindlassParameters.m)
if strcmp(subject,'s1_Poggensee')
    H0 = 0.0724; % unloaded foot arch height (m)
    a = 0.1126; % length from calcn origin to mtj (m)
    b = 0.1047; % length from mtj to mtpj (m)
    phi0 = 1.6798; % angle between a and b (rad)
elseif strcmp(subject,'subject1')
    H0 = 0.0724; % unloaded foot arch height (m)
    a = 0.1126;
    b = 0.1047;
    phi0 = 1.6798;
end

a = 0.08207;
b = 0.089638;
phi0 = 2.493499;
H0 = 0.027280;

L0 = sqrt(a^2 + b^2 - 2*a*b*cos(phi0));

% Get values resulting from Windlass mechanism, assuming infinitely stiff
% PF and no influence from other soft tissue.
l_0 = (1-cWL*(q_mtp*180/pi)/20)*L0; % foot arch length
h_0 = a*b./l_0*sin(phi0);
q_mt_0 = acos( (a^2 + b^2 - l_0.^2)/(2*a*b) ) - phi0;

% Get foot arch dimensions from mtj angle
phi = phi0 + q_mt;
l = sqrt(a^2 + b^2 - 2*a*b*cos(phi));
h = a*b./l.*sin(phi);

% Get linear PF stiffness
% Simplification, but order of magnitude of k is similar to
% doi:10.1016/j.clinbiomech.2004.06.002
% q1 = (-5:1:5)*pi/180; % use small range for tmt angle
% l1 = L0*(1 + cWLl*(q1+cos(q1)-1) );
% h1 = cWLh*q1 + H0;
% k_PF = nanmean(kMT_PF*q1./(h1.*(l1-L0)));

k_PF = 7.1994e+05*kMT_PF/1000;

dl_0 = 2e-3;

% Calculate moments
dl_PF = (l-l_0) .*( tanh( (l-l_0)*1e6 )+1 )/2; % PF elongation (>=0)
% F_PF = k_PF*dl_PF;
F_PF = k_PF*(dl_PF - dl_0*tanh(dl_PF/dl_0));
M_PF = F_PF*h; % moment from PF elongation

M_li = kTMT_li*q_mt; % moment from ligament elongation
M_d = dMT*qdot_mt; % moment from viscous friction

% Total passive moment
M = -(M_PF + M_li + M_d);

% Return some more outputs if the function is called from postprocessing
if nargout > 1
    varargout = {M_PF,F_PF,M_li,M_d,l,l_0,L0,h,h_0,H0,q_mt_0};
end

end