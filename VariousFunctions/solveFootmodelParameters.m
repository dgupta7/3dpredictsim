% Find the parameters for a new foot model, such that it is equivalent to
% the old footmodel, and to a more detailed footmodel.
% Every joint is taken at its default angle of 0°. This means that the 0°
% angle of the newly introduced tarsometarsal joint corresponds to its
% rigid position in the old model.

% note: These functions are not needed to run any simulation, their only
% purpose is to calculate parameters needed for the source code of the
% external functions.

clear all
close all
clc

%% parameters from old model
m_cmf_0 = 0.9688202167589921;
I_cmf_0 = [0.000906321, 0.00252475, 0.00265422]';
com_cmf_0 = [0.0913924, 0.0274177, 0]';
% ankle0 = [0 -0.41085882914747662 0]';  % in tibia ref
% subtalar0 = [-0.044572100000000003 -0.038339100000000001 -0.0072382799999999997]';  % in talus ref
mtp_0 = [0.163409678774199 -0.00182784875586352 0.000987038328166303]';   % in calcn ref
% derived params
mtp_02com_0 = com_cmf_0 - mtp_0;

%% parameters detailed model
com_calcn = [-0.0183999 -0.0127205 0.00168919]';
com_midfoot = [0.0136264 -0.00541573 -0.00226789]';
com_forefoot = [0.030362 -0.00692578 -0.000413788]';
chopart = [0.0221509217480428 0.00563051896570316 -0.000938103987425259]';   % in calcn ref
tarsometarsal = [0.0224650598180872 -0.0131925031448367 0.00457222612063733]';   % in midfoot ref
mtp = [0.0604862357021971 -0.0140878323268576 0.00286827055184947]';     % in forefoot ref
m_calcn = 0.289923030444027;
m_midfoot = 0.139888239487345;
m_forefoot = 0.236187108315951;
% derived params
tarsometarsal2com_calcn = -(chopart + tarsometarsal - com_calcn);
tarsometarsal2com_midfoot = -(tarsometarsal - com_midfoot);
tarsometarsal2com_cm = (tarsometarsal2com_calcn * m_calcn + tarsometarsal2com_midfoot * m_midfoot)/(m_calcn + m_midfoot);
mtp2com_cm = tarsometarsal2com_cm - mtp;
mtp2com_forefoot = com_forefoot - mtp;
mtp2com_cmf = (mtp2com_cm * (m_calcn + m_midfoot) + mtp2com_forefoot * m_forefoot)/(m_calcn + m_midfoot + m_forefoot);

% 3 relations about inertias needed: Open the forefoot of the detailed
% footmodel in Fusion 360. Convert mesh to solid, and calculate its
% inertias and mass for a unit mass density (assume constant mass distribution)

m_ff = 45.16;    % [g]
I_ff = [2.351, 4.183, 2.718]*1e4;  % [g mm^2]

I_ff0 = I_ff/m_ff *1e-6;  % [kg m^2 /kg] mass invariant inertia


% Scale factor to account for the fact that the detailed model is based on
% a different person.
sf = mtp_02com_0./mtp2com_cmf;

% initial guess (c = calcaneus, f = forefoot)
mc = m_cmf_0*0.7;
mf = m_cmf_0*0.3;
Icx = I_cmf_0(1)*0.7;
Icy = I_cmf_0(2)*0.7;
Icz = I_cmf_0(3)*0.7;
Ifx = I_cmf_0(1)*0.2;
Ify = I_cmf_0(2)*0.2;
Ifz = I_cmf_0(3)*0.2;
COMcx = com_cmf_0(1)*0.5;
COMcy = com_cmf_0(2)*0.9;
COMcz = com_cmf_0(3);
COMfx = com_cmf_0(1)*1.5;
COMfy = com_cmf_0(2)*1.1;
COMfz = com_cmf_0(3);
MTJx = com_cmf_0(1);
MTJy = com_cmf_0(2);
MTJz = com_cmf_0(3);

init_vars = [mc,mf,Icx,Icy,Icz,Ifx,Ify,Ifz,COMcx,COMcy,COMcz,COMfx,COMfy,COMfz,MTJx,MTJy,MTJz]';

% solve

optim_options = optimset('Display','off');

[x, fval, exitflag]=fsolve('footmodel',init_vars,optim_options,... 
    m_cmf_0,I_cmf_0,com_cmf_0,mtp_0,m_calcn,m_midfoot,m_forefoot,tarsometarsal2com_cm,com_forefoot,I_ff0,sf);
    
    if (exitflag ~= 1)
        disp 'The fsolve exit flag was not 1, probably no convergence!'
        exitflag
    end
    
% results for right foot
mc = x(1);
mf = x(2);
Ic = x(3:5);
If = x(6:8);
COMc = x(9:11);
COMf = x(12:14);
MTJ = x(15:17);
MTPJ = mtp_0-MTJ;


disp('I0    Ic    If')
[I_cmf_0 Ic If]

%% plot results in calcaneus ref frame
% This is the ref frame attached to the calcn in both the old model and the
% new model. The detailed model needs to be scaled.

% transform new model to ref frame
c2COMf = MTJ+COMf;
c2MTPJ = MTJ+MTPJ; %=mtp0
% transform detailed model to ref frame
c2mtpj = mtp_0;
c2mtj = mtp_0 - mtp;
c2com_ff = mtp_0 + mtp2com_forefoot;
c2com_cm = mtp_0 + mtp2com_cm;
% scale detailed model
c2mtpj_s = mtp_0; 
c2mtj_s = mtp_0 - mtp.*sf;
c2com_ff_s = mtp_0 + mtp2com_forefoot.*sf;
c2com_cm_s = mtp_0 + mtp2com_cm.*sf;


figure
subplot(2,1,1)
plot(mtp_0(1),mtp_0(2),'or')
hold on
grid on
plot(MTJ(1),MTJ(2),'or')
plot(c2COMf(1),c2COMf(2),'xr')
plot(COMc(1),COMc(2),'xr')
plot(com_cmf_0(1),com_cmf_0(2),'xg')
plot(c2mtpj(1),c2mtpj(2),'ob')
plot(c2mtj(1),c2mtj(2),'ob')
plot(c2com_ff(1),c2com_ff(2),'xb')
plot(c2com_cm(1),c2com_cm(2),'xb')
plot(c2mtpj_s(1),c2mtpj_s(2),'ok')
plot(c2mtj_s(1),c2mtj_s(2),'ok')
plot(c2com_ff_s(1),c2com_ff_s(2),'xk')
plot(c2com_cm_s(1),c2com_cm_s(2),'xk')
plot(0,0,'.')
title('side view')
subplot(2,1,2)
plot(mtp_0(1),mtp_0(3),'or')
hold on
grid on
plot(MTJ(1),MTJ(3),'or')
plot(c2COMf(1),c2COMf(3),'xr')
plot(COMc(1),COMc(3),'xr')
plot(com_cmf_0(1),com_cmf_0(3),'xg')
plot(c2mtpj(1),c2mtpj(3),'ob')
plot(c2mtj(1),c2mtj(3),'ob')
plot(c2com_ff(1),c2com_ff(3),'xb')
plot(c2com_cm(1),c2com_cm(3),'xb')
plot(c2mtpj_s(1),c2mtpj_s(3),'ok')
plot(c2mtj_s(1),c2mtj_s(3),'ok')
plot(c2com_ff_s(1),c2com_ff_s(3),'xk')
plot(c2com_cm_s(1),c2com_cm_s(3),'xk')
plot(0,0,'.')
title('top view')






