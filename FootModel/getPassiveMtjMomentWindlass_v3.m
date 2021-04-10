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
% R_mtth = 7.5e-3; % average radius of the metatarsal head
R_mtth = 9.5e-3; % radius of the first metatarsal head
l_toe = 0; % distance from metatarsal head to PF attachment point at toe

sf_PF = 1;
sf_li = 1;

%% Get specific parameters
if isempty(varargin)
    dMT = 0;
    nl = 1; % use nonlinear ligament stiffness
    kMT_li = 90; % linear stiffness
else
    S = varargin{1};
    dMT = S.dMT;
    nl = S.MT_li_nonl;
    kMT_li = S.kMT_li;
    if isfield(S,'sf_PF') && ~isempty(S.sf_PF)
        sf_PF = S.sf_PF;
    end
    if isfield(S,'sf_li') && ~isempty(S.sf_li)
        sf_li = S.sf_li;
    end
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
F_PF = f_PF_stiffness(l_PF)*sf_PF;
if nargout > 1
    % convert SX to double when called from post-processing
    F_PF = full(F_PF);
end

M_PF = -F_PF*MA_PF;

% other elastic structures
if nl
%     dl_0 = 3*pi/180;
%     dq_0 = 2*pi/180;
%     M_li = -kMT_li*( (q_mt-dq_0) - dl_0*tanh((q_mt-dq_0)/dl_0/1.2)) +... % toe-in and linear part
%             (-exp(20*(q_mt-20*pi/180)) + exp(-25*(q_mt+10*pi/180)))/2; % stiffening at the end
    
M_li = -9*(exp(4*(q_mt-2*pi/180))-1)*1.2 + 2*exp(-10*(q_mt+0.1));

%     M_li = -kMT_li*q_mt*(exp(5*(q_mt-10*pi/180)) + exp(-3*(q_mt+20*pi/180)))/2; % approx 1
    

%     k1 = 12;
%     t1 = 8*pi/180*sf_li;
%     f2 = 1.5*(2-sf_li);
%     k2 = 5;
%     t2 = 10*pi/180*sf_li;
%     M_li = (-exp(k1*(q_mt-t1)) + f2*exp(-k2*(q_mt+t2)))*kMT_li*5/90; % approx 2
%     M_li = (-exp(k1*(q_mt-t1)) + f2*exp(-k2*(q_mt+t2)))*5*3-2; % approx 2a
    
%     sf = 0.1;
%     t = 0.5;
%     c1 = 12*(2-t)/(1-t+sf);
%     c2 = 8*sf;
%     c3 = 1.5;
%     c4 = 5;
%     c5 = 10*sf;
%     M_li = (-exp(c1*(q_mt-c2*pi/180)) + c3*exp(-c4*(q_mt+c5*pi/180)))*kMT_li*5/90*(1+sf)/2;

%     sf = 0.3;
%     c1 = 12/(sf+0.5)*1.5;
%     c2 = 8*(1+sf)/2;
%     c3 = 1.5*(sf);%+0.4)/1.4;
%     c4 = 5;%/(sf+0.5)*1.5;
%     c5 = 10;%*(1+sf)/2;

% sf = 0.1;
% c1 = 12/(sf+1)*2;
% c2 = 8*(2+sf)/3;
% c3 = 1.5*(sf+0.5)/1.5;
% c4 = 5;%/(sf+5)*6;
% c5 = 10;%*(1+sf)/2;

% c1 = 25;
% c2 = 1;
% c3 = 3;
% c4 = 5;
% c5 = 1;
% sf=1;
% 
%     M_li = (-exp(c1*(q_mt-c2*pi/180)) + c3*exp(-c4*(q_mt+c5*pi/180)))*5/sf;

% M_li = (-exp(12*(q_mt*3-8*pi/180)) + 2*exp(-5*(q_mt*2+10*pi/180)))*5*2;

% M_li = (-exp(12*(q_mt*2-8*pi/180)) + 1.5*exp(-5*(q_mt*1.5+10*pi/180)))*5*5; % v5

% M_li = 0;

% M_li = -exp(25*(q_mt-5*pi/180)) + 2*exp(-15*(q_mt+5*pi/180)) +0.5;

% M_li = -2*exp(10*(q_mt-5*pi/180)) + 2*exp(-15*(q_mt+5*pi/180));

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