function [res] = StaticStanding(vars,S,F,f_AllPassiveTorques,Qmtp)
%     %% Load external functions
%     import casadi.*
%     % The external function performs inverse dynamics through the
%     % OpenSim/Simbody C++ API. This external function is compiled as a dll from
%     % which we create a Function instance using CasADi in MATLAB. More details
%     % about the external function can be found in the documentation.
%     pathmain        = pwd;
%     [pathRepo,~,~]  = fileparts(pathmain);
%     addpath(genpath(pathRepo));
%     % Loading external functions.
%     pathExternalFunctions = [pathRepo,'/ExternalFunctions'];
%     cd(pathExternalFunctions)
%     F  = external('F',S.ExternalFunc);
%     cd(pathmain);

    %% Indices external function
    % Indices of the elements in the external functions
    % External function: F
    % First, joint torques.
    jointi = getJointi_tmt();

    % Vectors of indices for later use
    residualsi          = jointi.pelvis.tilt:jointi.elb.r; % all
    % Number of degrees of freedom for later use
    nq.all      = length(residualsi); % all

%     %% CasADi functions
%     % We create several CasADi functions for later use
%     pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
%     PathDefaultFunc = fullfile(pathCasADiFunctions,S.CasadiFunc_Folders);
%     f_AllPassiveTorques = Function.load(fullfile(PathDefaultFunc,'f_AllPassiveTorques'));

    %%
    Qs_nsc = zeros(nq.all,1);
    Qdots_nsc = zeros(nq.all,1);
    A_nsc = zeros(nq.all,1);

    Qs_nsc(jointi.pelvis.ty) = vars(1);
    Qs_nsc(jointi.ankle.l) = vars(2);
    Qs_nsc(jointi.ankle.r) = vars(3);
    Qs_nsc(jointi.subt.l) = vars(4);
    Qs_nsc(jointi.subt.r) = vars(5);
    Qs_nsc(jointi.tmt.l) = vars(6);
    Qs_nsc(jointi.tmt.r) = vars(7);
    
    Qs_nsc(jointi.mtp.l) = Qmtp*pi/180;
    Qs_nsc(jointi.mtp.r) = Qmtp*pi/180;

    QsQdots_nsc = zeros(nq.all*2,1);
    QsQdots_nsc(1:2:end,:) = Qs_nsc;
    QsQdots_nsc(2:2:end,:) = Qdots_nsc;

    %%
    % Get passive joint torques
    Tau_passj_all = full(f_AllPassiveTorques(Qs_nsc(:,1),Qdots_nsc(:,1)));
    Tau_passj.ankle.l = Tau_passj_all(9);
    Tau_passj.ankle.r = Tau_passj_all(10);
    Tau_passj.subt.l = Tau_passj_all(11);
    Tau_passj.subt.r = Tau_passj_all(12);
    Tau_passj.tmt.l = Tau_passj_all(13);
    Tau_passj.tmt.r = Tau_passj_all(14);

    %%
    [Tj] = full(F([QsQdots_nsc(:,1);A_nsc(:,1)]));
    % Pelvis around resting height
    res(1) = tanh((Qs_nsc(jointi.pelvis.ty)-1.2*S.IG_PelvisY)*100)...
        -tanh((Qs_nsc(jointi.pelvis.ty)-0.8*S.IG_PelvisY)*100)+2;
    % Ankle, left
    res(2) = Tj(jointi.ankle.l,1)-(Tau_passj.ankle.l);
    % Ankle, right
    res(3) = Tj(jointi.ankle.r,1)-(Tau_passj.ankle.r);
    % Subtalar, left
    res(4) = Tj(jointi.subt.l,1)-(Tau_passj.subt.l);
    % Subtalar, right
    res(5) = Tj(jointi.subt.r,1)-(Tau_passj.subt.r );
    % Tmt, left
    res(6) = Tj(jointi.tmt.l,1) - Tau_passj.tmt.l;
    % Tmt, right
    res(7) = Tj(jointi.tmt.r,1) - Tau_passj.tmt.r;

end