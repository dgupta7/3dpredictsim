function [M, varargout] = getPassiveTmtjMomentWindlass(q_tmt,qdot_tmt,q_mtp,kTMT_li,kTMT_PF,dTMT,subject,cWL)
% This function returns the passive moment around the tarsometatarsal joint
% q_tmt: angle of tmt joint (rad)
% q_mtp: angle of mtp joint(rad)
% qdot_mtp: angular velocity of tmt joint (rad/s)
% kTMT_li: angular stiffness of tmt joint, from ligaments (Nm/rad)
% kTMT_PA: angular stiffness of tmt joint, from plantar fascia (Nm/rad)
%           Formulated like this for consistency with earlier models
% dTMT: angular viscous friction of tmt joint (Nms/rad)
% Subject: test subject
% cWL: relation between mtp angle and foot arch length
%           (see doi.org/10.1098/rsif.2018.0270)


% Get subject-specific constants, derived from foot geometry 
% (see \VariousFunctions\WindlassParameters.m)
if strcmp(subject,'s1_Poggensee')
    cWLq = -12.3342;
    cWLl = 0.2293;
    cWLh = -0.0364;
    L0 = 0.1628; % unloaded foot arch length (m)
    H0 = 0.0374; % unloaded foot arch height (m)
else
    
end

% Get values resulting from Windlass mechanism, assuming infinitely stiff
% PF and no influence from other soft tissue.

l_0 = (1-cWL*(q_mtp*180/pi)/20)*L0; % foot arch length
q_tmt_0 = cWL*cWLq*q_mtp; % tmt angle
h_0 = cWLh*q_tmt_0 + H0;

% Get foot arch dimensions from tmt angle
l = L0*(1 + cWLl*(q_tmt+cos(q_tmt)-1) );
h = cWLh*q_tmt + H0;

% Get linear PF stiffness
% Simplification, but order of magnitude of k is similar to
% doi:10.1016/j.clinbiomech.2004.06.002
q1 = (-5:1:5)*pi/180; % use small range for tmt angle
l1 = L0*(1 + cWLl*(q1+cos(q1)-1) );
h1 = cWLh*q1 + H0;
k_PF = nanmean(kTMT_PF*q1./(h1.*(l1-L0)));

% Calculate moments
dl_PF = (l-l_0) .*( tanh( (l-l_0)*1e6 )+1 )/2; % PF elongation (>=0)
F_PF = k_PF*dl_PF;
M_PF = F_PF*h; % moment from PF elongation

M_li = kTMT_li*q_tmt; % moment from ligament elongation
M_d = dTMT*qdot_tmt; % moment from viscous friction

% Total passive moment
M = -(M_PF + M_li + M_d);

% Return some more outputs if the function is called from postprocessing
if nargout > 1
    varargout = {M_PF,F_PF,M_li,M_d,l,l_0,L0,h,h_0,H0,q_tmt_0};
end


end