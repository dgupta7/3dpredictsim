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
R_mtth = 7.5e-3; % average radius of the metatarsal head
% R_mtth = 9.5e-3; % radius of the first metatarsal head

sf_PF = 1; % scale factor for PF stiffness
mtj_stiffness = 'Ker1987'; % model for resulting foot stiffness
dMT = 0;
nl = 1; % use nonlinear ligament stiffness
kMT_li = 90; % linear stiffness

%% Get specific parameters
if isempty(varargin)
    
else
    S = varargin{1};
    if isfield(S,'dMT') && ~isempty(S.dMT)
        dMT = S.dMT;
    end
    if isfield(S,'MT_li_nonl') && ~isempty(S.MT_li_nonl)
        nl = S.MT_li_nonl;
    end
    if isfield(S,'kMT_li') && ~isempty(S.kMT_li)
        kMT_li = S.kMT_li;
    end
    if isfield(S,'sf_PF') && ~isempty(S.sf_PF)
        sf_PF = S.sf_PF;
    end
    if isfield(S,'mtj_stiffness') && ~isempty(S.mtj_stiffness)
        mtj_stiffness = S.mtj_stiffness;
    end
end


%% Geometry Windlass mechanism
% based on midtarsal joint angle
phi = phi0 + q_mt; % top angle of WL triangle
l_PF_fa = sqrt(calcnPF2mtj^2 + mtj2mttPF^2 - 2*calcnPF2mtj*mtj2mttPF*cos(phi)); % length of PF spanning arch
MA_PF = calcnPF2mtj*mtj2mttPF/l_PF_fa*sin(phi); % moment arm of PF to mtj

% Plantar fascia length
l_PF = l_PF_fa + R_mtth*(pi/2+q_mtp) -0.0017; % constant term to correct path length to physical length

%% Torques
% plantar fascia
try
    F_PF = f_PF_stiffness(l_PF)*sf_PF;
    if nargout > 2
        % convert SX to double when called from post-processing
        F_PF = full(F_PF);
    end
catch
%     warning('No valid function to describe the plantar fascia stiffness model, using 0 instead.')
    F_PF = 0;
end

M_PF = -F_PF*MA_PF;

% other elastic structures
if nl
    if strcmp(mtj_stiffness,'Gefen2001')
        M_li = -9*(exp(4*(q_mt-2*pi/180))-1)*1.2 + 2*exp(-10*(q_mt+0.1));
        
    elseif strcmp(mtj_stiffness,'Ker1987')
        % calculated in ligaments_torques_Ker87_v2.m
        M_li = 2.16362 + -55.325857*q_mt^1 + -614.60128*q_mt^3 + 1732.1067*q_mt^5 + ...
            7091.9679*q_mt^7 + -24459.968*q_mt^9 -1;
            
    elseif strcmp(mtj_stiffness,'fitted')
%         % for Natali2010, with ls=144
%         t1 = 0*pi/180;
%         c1 = 10;
%         c2 = 20;
% 
%         t2 = 10*pi/180;
%         c3 = 1;
%         c4 = 50;
%         c5 = 28;

        t1 = 3*pi/180;
        c1 = 10;
        c2 = 25;

        t2 = 10*pi/180;
        c3 = 2;
        c4 = 40;
        c5 = 8;

        M_li = -c1*exp(c2*(q_mt-t1)) + c3*exp(-c4*(q_mt+t2)) + c5;
        
    end
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