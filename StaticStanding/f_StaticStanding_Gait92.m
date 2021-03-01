function [] = f_StaticStanding_Gait92(S)

AddCasadiPaths();

%% Default settings
S = GetDefaultSettings(S);

%% Load external functions
import casadi.*
% The external function performs inverse dynamics through the
% OpenSim/Simbody C++ API. This external function is compiled as a dll from
% which we create a Function instance using CasADi in MATLAB. More details
% about the external function can be found in the documentation.
pathmain        = pwd;
[pathRepo,~,~]  = fileparts(pathmain);
addpath(genpath(pathRepo));
% Loading external functions.
pathExternalFunctions = [pathRepo,'/ExternalFunctions'];
cd(pathExternalFunctions)
F  = external('F',S.ExternalFunc);
F1 = external('F',S.ExternalFunc2);
cd(pathmain);

%% CasADi functions
% We create several CasADi functions for later use
pathCasADiFunctions = [pathRepo,'/CasADiFunctions'];
PathDefaultFunc = fullfile(pathCasADiFunctions,S.CasadiFunc_Folders);
f_AllPassiveTorques = Function.load(fullfile(PathDefaultFunc,'f_AllPassiveTorques'));
    
%% Indices external function
% Indices of the elements in the external functions
% External function: F
% First, joint torques.
jointi = getJointi_tmt();
% Vectors of indices for later use
residualsi          = jointi.pelvis.tilt:jointi.elb.r; % all
% Number of degrees of freedom for later use
nq.all      = length(residualsi); % all
% External function: F1
GRFi.r      = 34:36;
GRFi.l      = 37:39;

%%
% mtp angles to be considered
Qs_mtp = [-30:5:45];
n_mtp = length(Qs_mtp);

% initial guess for solver
vars_init = [S.IG_PelvisY; zeros(6,1)];
optim_options = optimset('Display','off');

% Windlass parameters
kTMT_li = 1.5/(pi/180)/5;
kTMT_PF = S.kTMT;
dTMT = S.dTMT;
cWL = S.cWL;

% Declare arrays for postprocessing results
Qs = zeros(nq.all,n_mtp);
Qdots = zeros(nq.all,n_mtp);
Qddots = zeros(nq.all,n_mtp);
QsQdots = zeros(nq.all*2,n_mtp);
GRF_r = zeros(3,n_mtp);
M = zeros(n_mtp,1);
M_PF = zeros(n_mtp,1);
F_PF = zeros(n_mtp,1);
l_fa = zeros(n_mtp,1);
h_fa = zeros(n_mtp,1);
l0_fa = zeros(n_mtp,1);
h0_fa = zeros(n_mtp,1);
q_tmt_0 = zeros(n_mtp,1);


for i=1:n_mtp
    %% solve the static situation
    [x,fval,~]=fsolve('StaticStanding',vars_init,optim_options,S,F,f_AllPassiveTorques,Qs_mtp(i));
    
    %% postprocess results
    % get states
    Qs(jointi.pelvis.ty,i) = x(1);
    Qs(jointi.ankle.l,i) = x(2);
    Qs(jointi.ankle.r,i) = x(3);
    Qs(jointi.subt.l,i) = x(4);
    Qs(jointi.subt.r,i) = x(5);
    Qs(jointi.tmt.l,i) = x(6);
    Qs(jointi.tmt.r,i) = x(7);
    Qs(jointi.mtp.l,i) = Qs_mtp(i)*pi/180;
    Qs(jointi.mtp.r,i) = Qs_mtp(i)*pi/180;

    QsQdots(1:2:end,i) = Qs(:,i);
    QsQdots(2:2:end,i) = Qdots(:,i);

    % call external function
    [res] = full(F1([QsQdots(:,i);Qddots(:,i)]));
    % get GRF
    GRF_r(:,i) = res(GRFi.r);

    [Mi, M_PFi,F_PFi,~,~,li,l0i,L0,hi,h0i,H0,q_tmt_0i] = ...
            getPassiveTmtjMomentWindlass(Qs(jointi.tmt.r,i),0,Qs(jointi.mtp.r,i),...
            kTMT_li,kTMT_PF,dTMT,S.subject,cWL);

    M(i) = Mi;
    M_PF(i) = M_PFi;
    l_fa(i) = li;
    h_fa(i) = hi;
    F_PF(i) = F_PFi;
    l0_fa(i) = l0i;
    h0_fa(i) = h0i;
    q_tmt_0(i) = q_tmt_0i;
    
    
    
end


% save results
R.S = S;
R.Qs_mtp_deg = Qs_mtp;
R.Qs_rad = Qs;
R.GRF_r = GRF_r;
R.M_WL = M;
R.M_PF = M_PF;
R.F_PF = F_PF;
R.l_fa = l_fa;
R.l0_fa = l0_fa;
R.L0 = L0;
R.h_fa = h_fa;
R.h0_fa = h0_fa;
R.H0 = H0;
R.Q_tmt_0_rad = q_tmt_0;

OutFolder = fullfile(pathRepo,'Results',S.ResultsFolder);
FilenameAnalysis = fullfile(OutFolder,[S.savename '_pp.mat']);
save(FilenameAnalysis,'R');


%% make plots (move to separate function)
figure
subplot(2,3,1)
hold on
plot(Qs_mtp,l0_fa./L0,'--')
plot(Qs_mtp,l_fa./L0,'-')
xlabel('mtp angle (°)')
ylabel('arch length (-)')

subplot(2,3,2)
hold on
plot(Qs_mtp,h0_fa./H0,'--')
plot(Qs_mtp,h_fa./H0,'-')
xlabel('mtp angle (°)')
ylabel('arch height (-)')

subplot(2,3,3)
hold on
plot(Qs_mtp,GRF_r(2,:))
xlabel('mtp angle (°)')
ylabel('GRF_y (N)')

subplot(2,3,4)
hold on
plot(Qs_mtp,F_PF)
xlabel('mtp angle (°)')
ylabel('F PF (N)')

subplot(2,3,5)
hold on
plot(Qs_mtp,q_tmt_0*180/pi,'--')
plot(Qs_mtp,Qs(jointi.tmt.r,:)*180/pi,'-')
xlabel('mtp angle (°)')
ylabel('tmt angle (°)')







end