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

sf_PF = 1; % scale factor for PF stiffness
mtj_stiffness = 'Ker1987'; % model for resulting foot stiffness
dMT = 0;
nl = 0; % use nonlinear ligament stiffness
kMT_li = 300; % linear stiffness
kMT_li2 = kMT_li;
k_sta = 0; % externally added stiffness to the mtj (e.g. shoe, insole)
poly = 1; % polynomial approximation instead of sin,cos

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
    if isfield(S,'kMT_li2') && ~isempty(S.kMT_li2)
        kMT_li2 = S.kMT_li2;
    else
        kMT_li2 = kMT_li;
    end
    if isfield(S,'sf_PF') && ~isempty(S.sf_PF)
        sf_PF = S.sf_PF;
    end
    if isfield(S,'mtj_stiffness') && ~isempty(S.mtj_stiffness)
        mtj_stiffness = S.mtj_stiffness;
    end
    if isfield(S,'stiffen_arch') && ~isempty(S.stiffen_arch)
        k_sta = S.stiffen_arch;
    end
    if isfield(S,'R_mtth') && ~isempty(S.R_mtth)
        R_mtth = S.R_mtth;
    end
    if isfield(S,'WLpoly') && ~isempty(S.WLpoly)
        poly = S.WLpoly;
    end
end
F_PIM = 0;
if length(varargin)>=2
    if isfield(S,'PIM') && ~isempty(S.PIM) && S.PIM==1
        F_PIM = varargin{2};
    else
        sf_PF = varargin{2}*20;
    end
end


%% Geometry Windlass mechanism
if poly == 0
    % based on midtarsal joint angle
    % top angle of WL triangle
    phi = phi0 + q_mt;
    % length of PF spanning arch
    l_PF_fa = sqrt(calcnPF2mtj^2 + mtj2mttPF^2 - 2*calcnPF2mtj*mtj2mttPF*cos(phi));
    % moment arm of PF to mtj
    MA_PF = calcnPF2mtj*mtj2mttPF/l_PF_fa*sin(phi);
    
    % Plantar fascia length
    l_PF = l_PF_fa + R_mtth*(pi/4+q_mtp) + 0.0042; % constant term to correct to physical length

else
%     % see \FootModel\windlassGeometryPolynomials.m
%     % length of PF spanning arch
    l_PF_fa = 0.1392179 + 0.0374482.*q_mt.^1 + -0.0166876.*q_mt.^2 + ...
        -0.001758651.*q_mt.^3 + 0.0004480769.*q_mt.^4;
%     % moment arm of PF to mtj
%     MA_PF = 0.0374478 + -0.03337403.*q_mt.^1 + -0.005255987.*q_mt.^2 + ...
%         0.001767266.*q_mt.^3 + -0.0001071423.*q_mt.^4 + 9.858065e-05.*q_mt.^5;

    [MA_PF,l_PF,~] = getPlantarFasciaLengthVelocity(q_mt,0,q_mtp,0,R_mtth);
    
end

    
%% Torques
% plantar fascia
try
    F_PF = f_PF_stiffness(l_PF);
    if nargout > 2
        % convert SX to double when called from post-processing
        F_PF = full(F_PF);
    end
catch
%     warning('No valid function to describe the plantar fascia stiffness model, using 0 instead.')
    F_PF = 0;
end
F_PF = F_PF*sf_PF;
F_PF_tot = F_PF + F_PIM;   
M_PF = -F_PF_tot*MA_PF;

% other elastic structures
if nl
    if strcmp(mtj_stiffness,'Gefen2002')
        M_li = -9*(exp(4*(q_mt-2*pi/180))-1)*1.2 + 2*exp(-10*(q_mt+0.1));
        
    elseif strcmp(mtj_stiffness,'Ker1987')
        % calculated in ligaments_torques_Ker87_v2.m
        M_li = 2.16362 + -55.325857*q_mt^1 + -614.60128*q_mt^3 + 1732.1067*q_mt^5 + ...
            7091.9679*q_mt^7 + -24459.968*q_mt^9 -1;
        
    elseif strcmp(mtj_stiffness,'Song2011')
        % S. Song, C. LaMontagna, S. H. Collins en H. Geyer,
        % „The Effect of Foot Compliance Encoded in the Windlass Mechanism 
        % on the Energetics of Human Walking,” in 35th Annual International 
        % Conference of the IEEE EMBS, Osaka, Japan, 2013. 
        M_li = -800*q_mt;
        
    elseif strcmp(mtj_stiffness,'signed_lin')
        y_p = -kMT_li*q_mt.*(tanh(q_mt*100-0.5)+1)/2;
        y_n = -kMT_li2*q_mt.*(-tanh(q_mt*100+0.5)+1)/2;
        M_li = y_p + y_n  + 1.55;
        
    elseif strcmp(mtj_stiffness,'fitted1')
        % for Natali2010, with ls=148
        t1 = 3*pi/180;
        c1 = 10;
        c2 = 25;

        t2 = 10*pi/180;
        c3 = 2;
        c4 = 40;
        c5 = 8;
        
        M_li = -c1*exp(c2*(q_mt-t1)) + c3*exp(-c4*(q_mt+t2)) + c5;
        
    elseif strcmp(mtj_stiffness,'fitted2')
        % v2
        t1 = 5*pi/180;
        c1 = 10;
        c2 = 25;

        t2 = 10*pi/180;
        c3 = 2;
        c4 = 40;
        c5 = 1;

        M_li = -c1*exp(c2*(q_mt-t1)) + c3*exp(-c4*(q_mt+t2)) + c5;
        
    elseif strcmp(mtj_stiffness,'fitted3')
        % v3 same as 1, but with influence from mtp
        t1 = 3*pi/180;
        c1 = 10;
        c2 = 25;

        t2 = 10*pi/180;
        c3 = 2;
        c4 = 40;
        c5 = 8;
        
        c6 = 10;
        
        M_li = -c1*exp(c2*(q_mt-t1)) + c3*exp(-c4*(q_mt+t2)) + c5 - c6*q_mtp;
        
        elseif strcmp(mtj_stiffness,'fitted4')
            t1 = 5*pi/180;
            c1 = 10;
            c2 = 25;

            t2 = 10*pi/180;
            c3 = 2;
            c4 = 40;
            c5 = 8;
            
            M_li = -c1*exp(c2*(q_mt-t1)) + c3*exp(-c4*(q_mt+t2)) + c5;
            
        elseif strcmp(mtj_stiffness,'fitted5')
            t1 = 5*pi/180;
            c1 = 10;
            c2 = 25;

            t2 = 10*pi/180;
            c3 = 10;
            c4 = 15;
            c5 = 10;
            
            M_li = -c1*exp(c2*(q_mt-t1)) + c3*exp(-c4*(q_mt+t2)) + c5;
            
        elseif strcmp(mtj_stiffness,'fitted6')
            t1 = 5*pi/180;
            c1 = 10;
            c2 = 25;

            t2 = 10*pi/180;
            c3 = 5;
            c4 = 20;
            c5 = 4;
            
            M_li = -c1*exp(c2*(q_mt-t1)) + c3*exp(-c4*(q_mt+t2)) + c5;
            
    end
else
    M_li = -kMT_li*q_mt + 1.55;
end
       
% additional external stiffness
if k_sta ~= 0
    M_li = M_li - k_sta*q_mt;
end

% viscous damping
M_d = -dMT*qdot_mt;

% Total passive moment
M_mtj = M_PF + M_li + M_d;

% Reaction torque on toes
M_mtpj = -R_mtth*F_PF_tot;

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