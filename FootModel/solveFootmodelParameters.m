% Find the parameters for a new foot model (right foot), such that it is equivalent to
% the old footmodel, and to a more detailed footmodel.
% Every joint is taken at its default angle of 0°. This means that the 0°
% angle of the newly introduced m2.m.tmtj joint corresponds to its
% rigid position in the old model.

% note: These functions are not needed to run any simulation, their only
% purpose is to calculate parameters needed for the source code of the
% external functions.

% Parameter name explenation:
% 1) Model used
%   m0 = original model (0 dof)
%   m2 = detailed reference model (2 dof)
%   m1a = original model with tmt joint added in (1 dof)
%   m1b = original model with midtarsal joint added in (1 dof)
% 2) Body or reference frame
%   gnd = ground reference frame
%   c = calcaneus
%   m = midfoot (cuboid, navicular and cuneiform bones)
%   f = forefoot (metatarsi)
%   combination: 2 or 3 of the above are locked together in their defined
%   neutral position.
% 3) Variable (expressed in frame of 2)
%   m = mass (kg)
%   I = inertias [Ixx,Iyy,Izz] (kg m^2)
%   COM = centre of mass (m)
%   mtj, tmtk, mtpj = position of joint (m)
%   mtpj2cmfCOM = vector from mtp joint to the combined COM of c, m and f

% Author: Lars D'Hondt (march 2021)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
clc

%%
add_tmtj = 1;
add_mtj = 1;

%% parameters from original model
% Pog s1
m0.cmf.m = 0.9688202167589921;
m0.cmf.I = [0.000906321, 0.00252475, 0.00265422]';
m0.cmf.COM = [0.0913924, 0.0274177, 0]';
m0.cmf.mtpj = [0.163409678774199 -0.00182784875586352 0.000987038328166303]';   % in calcn ref
m0.t.subt = [-0.044572, -0.038339, 0.00723828]';

% % Fal s1
% m0.cmf.m = 0.938544584985273;
% m0.cmf.I = [0.000877997854457612, 0.00244585116598906, 0.00257127943091158]';
% m0.cmf.COM = [0.0913924, 0.0274177, 0]';
% m0.cmf.mtpj = [0.163409678774199, -0.00182784875586352, 0.000987038328166303]';   % in calcn ref


% location of contact spheres
% on calcn:
locSphere_1_r=[-0.00042152, -0.01, -0.0049972];
locSphere_2_r=[0.06, -0.01, 0.020001];
locSphere_3_r=[0.165, -0.01, 0.021183];
locSphere_4_r=[0.165, -0.01, -0.01];
% on toes:
locSphere_5_r=[0.053154, -0.01, -0.0034173];
locSphere_6_r=[1.7381e-06, -0.01, 0.022294];

    
% derived params
m0.gnd.mtpj2cmfCOM = m0.cmf.COM - m0.cmf.mtpj;
m0.gnd.talus2mtpj = m0.t.subt + m0.cmf.mtpj;

%% parameters detailed model
m2.c.COM = [-0.0183999 -0.0127205 0.00168919]';
m2.m.COM = [0.0136264 -0.00541573 -0.00226789]';
m2.f.COM = [0.030362 -0.00692578 -0.000413788]';
m2.c.mtj = [0.0221509217480428 0.00563051896570316 -0.000938103987425259]';   % in calcn ref
m2.m.tmtj = [0.0224650598180872 -0.0131925031448367 0.00457222612063733]';   % in midfoot ref
m2.f.mtpj = [0.0604862357021971 -0.0140878323268576 0.00286827055184947]';     % in forefoot ref
m2.c.m = 0.289923030444027;
m2.m.m = 0.139888239487345;
m2.f.m = 0.236187108315951;
m2.t.subt = [0.010347 -0.03163 0.00134647]';
% derived params
m2.gnd.tmtj2cCOM = -m2.c.mtj - m2.m.tmtj + m2.c.COM;
m2.gnd.tmtj2mCOM = -m2.m.tmtj + m2.m.COM;
m2.gnd.tmtj2cmCOM = (m2.gnd.tmtj2cCOM * m2.c.m + m2.gnd.tmtj2mCOM * m2.m.m)/(m2.c.m + m2.m.m);
m2.gnd.mtj2cCOM = -m2.c.mtj + m2.c.COM;
m2.gnd.mtpj2cCOM = m2.gnd.tmtj2cCOM - m2.f.mtpj;
m2.gnd.mtpj2cmCOM = m2.gnd.tmtj2cmCOM - m2.f.mtpj;
m2.gnd.mtpj2fCOM = -m2.f.mtpj + m2.f.COM;
m2.gnd.mtpj2cmfCOM = (m2.gnd.mtpj2cmCOM * (m2.c.m + m2.m.m) + m2.gnd.mtpj2fCOM * m2.f.m)/(m2.c.m + m2.m.m + m2.f.m);
m2.gnd.mtpj2mfCOM = (m2.gnd.mtpj2cmfCOM * (m2.c.m + m2.m.m + m2.f.m) - m2.gnd.mtpj2cCOM * m2.c.m)/(m2.m.m + m2.f.m);
m2.gnd.mtj2mfCOM = m2.gnd.mtpj2mfCOM + m2.f.mtpj + m2.m.tmtj;
m2.gnd.talus2mtpj = m2.t.subt + m2.c.mtj + m2.m.tmtj + m2.f.mtpj;

% Scale factor to account for the fact that the detailed model is based on
% a different person.
% Distance from mtpj to com of calcn + midfoot + forefoot, since this vector can be specified in all 3 models.
sf = m0.gnd.mtpj2cmfCOM./m2.gnd.mtpj2cmfCOM;
% sf = m0.gnd.talus2mtpj./m2.gnd.talus2mtpj;

% transform detailed model to ref frame
m2.gnd.mtpj = m0.cmf.mtpj;
m2.gnd.tmtj = m0.cmf.mtpj - m2.f.mtpj;
m2.gnd.mtj = m0.cmf.mtpj - m2.f.mtpj - m2.m.tmtj;
m2.gnd.fCOM = m0.cmf.mtpj + m2.gnd.mtpj2fCOM;
m2.gnd.cmCOM = m0.cmf.mtpj + m2.gnd.mtpj2cmCOM;
m2.gnd.mfCOM = m0.cmf.mtpj + m2.gnd.mtpj2mfCOM;
m2.gnd.cCOM = m0.cmf.mtpj + m2.gnd.mtpj2cCOM;
% scale detailed model
m2.gnd_s.mtpj = m0.cmf.mtpj; 
m2.gnd_s.tmtj = m0.cmf.mtpj - m2.f.mtpj.*sf;
m2.gnd_s.mtj = m0.cmf.mtpj - (m2.f.mtpj + m2.m.tmtj).*sf;
m2.gnd_s.fCOM = m0.cmf.mtpj + m2.gnd.mtpj2fCOM.*sf;
m2.gnd_s.cmCOM = m0.cmf.mtpj + m2.gnd.mtpj2cmCOM.*sf;
m2.gnd_s.mfCOM = m0.cmf.mtpj + m2.gnd.mtpj2mfCOM.*sf;
m2.gnd_s.cCOM = m0.cmf.mtpj + m2.gnd.mtpj2cCOM.*sf;
    

if add_tmtj
    %% approximate moments of inertia
    % Not enough relations to describe inertias => approximate some values
    % Open the forefoot of the detailed
    % footmodel in Fusion 360. Convert mesh to solid, specify a mass density that satisfies the total mass, and calculate its
    % moments of inertia at center of mass. Cross-elements are neglected.

    Ifx = 0.1893e-3;     %kg m^2 /kg
    Ify = 0.6381e-3;
    Ifz = 0.5466e-3;
    m1a.f.I = [Ifx; Ify; Ifz]*m0.cmf.m;


    %% solve
    % initial guess
    mcm = m0.cmf.m*0.7;
    mf = m0.cmf.m*0.3;
    Icmx = m0.cmf.I(1)*0.7;
    Icmy = m0.cmf.I(2)*0.7;
    Icmz = m0.cmf.I(3)*0.7;
    COMcmx = m0.cmf.COM(1)*0.5;
    COMcmy = m0.cmf.COM(2)*0.9;
    COMcmz = m0.cmf.COM(3);
    COMfx = m0.cmf.COM(1)*1.5;
    COMfy = m0.cmf.COM(2)*1.1;
    COMfz = m0.cmf.COM(3);
    TMTJx = m0.cmf.COM(1);
    TMTJy = m0.cmf.COM(2);
    TMTJz = m0.cmf.COM(3);

    init_vars = [mcm,mf,COMcmx,COMcmy,COMcmz,COMfx,COMfy,COMfz,TMTJx,TMTJy,TMTJz,Icmx,Icmy,Icmz]';

    % solver
    optim_options = optimset('Display','off');

    [x, fval,~]=fsolve('f_footmodel_add_tmtj',init_vars,optim_options,... 
        m0.cmf.m,m0.cmf.I,m0.cmf.COM,m2.c.m,m2.m.m,m2.f.m,m2.gnd.tmtj2cmCOM,m2.f.COM,m1a.f.I,sf);

    %% results for right foot
    % read results
    m1a.cm.m = x(1);
    m1a.f.m = x(2);
    m1a.cm.COM = x(3:5);
    m1a.f.COM = x(6:8);
    m1a.cm.tmtj = x(9:11);
    m1a.cm.I = x(12:14);
    
    % build more vectorrs
    m1a.f.mtpj = m0.cmf.mtpj-m1a.cm.tmtj;
    m1a.gnd.fCOM = m1a.cm.tmtj + m1a.f.COM;
    m1a.gnd.mtpj = m1a.cm.tmtj + m1a.f.mtpj;

    locSphere_3_r_new = locSphere_3_r - m1a.cm.tmtj';
    locSphere_4_r_new = locSphere_4_r - m1a.cm.tmtj';

    %% plot results in calcaneus ref frame
    % This is the ref frame attached to the calcn in both the old model and the
    % new model. The detailed model needs to be scaled.

    figure
    subplot(2,1,1)
    hold on
    grid on

    plot(m2.gnd.mtpj(1),m2.gnd.mtpj(2),'.b')
    plot(m2.gnd.tmtj(1),m2.gnd.tmtj(2),'.b')
    plot(m2.gnd.fCOM(1),m2.gnd.fCOM(2),'*b')
    plot(m2.gnd.cmCOM(1),m2.gnd.cmCOM(2),'*b')

    plot(m2.gnd_s.mtpj(1),m2.gnd_s.mtpj(2),'.k')
    plot(m2.gnd_s.tmtj(1),m2.gnd_s.tmtj(2),'.k')
    plot(m2.gnd_s.fCOM(1),m2.gnd_s.fCOM(2),'*k')
    plot(m2.gnd_s.cmCOM(1),m2.gnd_s.cmCOM(2),'*k')

    plot(m1a.gnd.mtpj(1),m1a.gnd.mtpj(2),'or')
    plot(m1a.cm.tmtj(1),m1a.cm.tmtj(2),'or')
    plot(m1a.gnd.fCOM(1),m1a.gnd.fCOM(2),'xr')
    plot(m1a.cm.COM(1),m1a.cm.COM(2),'xr')

    plot(m0.cmf.mtpj(1),m0.cmf.mtpj(2),'og')
    plot(m0.cmf.COM(1),m0.cmf.COM(2),'xg')
    plot(0,0,'og')
    title('side view (sagittal plane)')

    plot(locSphere_1_r(1),locSphere_1_r(2),'*c')
    plot(locSphere_2_r(1),locSphere_2_r(2),'*c')
    plot(locSphere_3_r(1),locSphere_3_r(2),'*c')
    plot(locSphere_4_r(1),locSphere_4_r(2),'*c')
    plot(m0.cmf.mtpj(1)+locSphere_5_r(1),m0.cmf.mtpj(2)+locSphere_5_r(2),'*c')
    plot(m0.cmf.mtpj(1)+locSphere_6_r(1),m0.cmf.mtpj(2)+locSphere_6_r(2),'*c')

    subplot(2,1,2)
    hold on
    grid on

    plot(m2.gnd.mtpj(1),m2.gnd.mtpj(3),'.b')
    plot(m2.gnd.tmtj(1),m2.gnd.tmtj(3),'.b')
    plot(m2.gnd.fCOM(1),m2.gnd.fCOM(3),'*b')
    plot(m2.gnd.cmCOM(1),m2.gnd.cmCOM(3),'*b')

    plot(m2.gnd_s.mtpj(1),m2.gnd_s.mtpj(3),'.k')
    plot(m2.gnd_s.tmtj(1),m2.gnd_s.tmtj(3),'.k')
    plot(m2.gnd_s.fCOM(1),m2.gnd_s.fCOM(3),'*k')
    plot(m2.gnd_s.cmCOM(1),m2.gnd_s.cmCOM(3),'*k')

    plot(m1a.gnd.mtpj(1),m1a.gnd.mtpj(3),'or')
    plot(m1a.cm.tmtj(1),m1a.cm.tmtj(3),'or')
    plot(m1a.gnd.fCOM(1),m1a.gnd.fCOM(3),'xr')
    plot(m1a.cm.COM(1),m1a.cm.COM(3),'xr')

    plot(m0.cmf.mtpj(1),m0.cmf.mtpj(3),'og')
    plot(m0.cmf.COM(1),m0.cmf.COM(3),'xg')
    plot(0,0,'og')
    title('bottom view')

    plot(locSphere_1_r(1),locSphere_1_r(3),'*c')
    plot(locSphere_2_r(1),locSphere_2_r(3),'*c')
    plot(locSphere_3_r(1),locSphere_3_r(3),'*c')
    plot(locSphere_4_r(1),locSphere_4_r(3),'*c')
    plot(m0.cmf.mtpj(1)+locSphere_5_r(1),m0.cmf.mtpj(3)+locSphere_5_r(3),'*c')
    plot(m0.cmf.mtpj(1)+locSphere_6_r(1),m0.cmf.mtpj(3)+locSphere_6_r(3),'*c')
    

    %% relating vectors to foot arch compression
    % based on DOI: 10.1038/srep19403

    a = m1a.cm.tmtj(1:2);
    b = m1a.f.mtpj(1:2);

    l_0 = norm(a+b);
    h0 = -b(2); 
    % Difference in y-coordinate between calcn origin and m1a.f.mtpj 
    % results in 1° difference, so it is omitted.

    c0 = acos(h0/norm(a));
    d0 = acos(h0/norm(b));

    tmt0 = (c0+d0)*180/pi;

    h1 = h0*0.8; 
    c1 = acos(h1/norm(a));
    d1 = acos(h1/norm(b));

    tmt1 = (c1+d1)*180/pi;

    tmt_bound = tmt1 - tmt0;
    % So 15° bound is sensible

    h2 = h0*0.87; 
    c2 = acos(h2/norm(a));
    d2 = acos(h2/norm(b));

    tmt2 = (c2+d2)*180/pi;

    tmt_2 = tmt2 - tmt0;
    
end
if add_mtj
    %% approximate moments of inertia
    % Not enough relations to describe inertias => approximate some values
    % Open the forefoot of the detailed
    % footmodel in Fusion 360. Convert mesh to solid, specify a mass density that satisfies the total mass, and calculate its
    % moments of inertia at center of mass. Cross-elements are neglected.

    Icx = 0.08842e-3;     %kg m^2 /kg (normalised with cmf mass)
    Icy = 0.1807e-3;
    Icz = 0.2077e-3;
    m1b.c.I = [Icx; Icy; Icz]*m0.cmf.m;


    %% solve
    % initial guess
    mc = m0.cmf.m*0.7;
    mmf = m0.cmf.m*0.3;
    Imfx = m0.cmf.I(1)*0.7;
    Imfy = m0.cmf.I(2)*0.7;
    Imfz = m0.cmf.I(3)*0.7;
    COMcx = m0.cmf.COM(1)*0.5;
    COMcy = m0.cmf.COM(2)*0.9;
    COMcz = m0.cmf.COM(3);
    COMmfx = m0.cmf.COM(1)*1.5;
    COMmfy = m0.cmf.COM(2)*1.1;
    COMmfz = m0.cmf.COM(3);
    MTJx = m0.cmf.COM(1);
    MTJy = m0.cmf.COM(2);
    MTJz = m0.cmf.COM(3);

    init_vars = [mc,mmf,COMcx,COMcy,COMcz,COMmfx,COMmfy,COMmfz,MTJx,MTJy,MTJz,Imfx,Imfy,Imfz]';

    % solver
    optim_options = optimset('Display','off');

    [x, fval,~]=fsolve('f_footmodel_add_mtj',init_vars,optim_options,... 
        m0.cmf.m,m0.cmf.I,m0.cmf.COM,m2.c.m,m2.m.m,m2.f.m,m2.gnd.mtj2mfCOM,m2.gnd.mtj2cCOM,m1b.c.I,sf);

    %% results for right foot
    % read results
    m1b.c.m = x(1);
    m1b.mf.m = x(2);
    m1b.c.COM = x(3:5);
    m1b.mf.COM = x(6:8);
    m1b.c.mtj = x(9:11);
    m1b.mf.I = x(12:14);
    
    % build more vectorrs
    m1b.mf.mtpj = m0.cmf.mtpj - m1b.c.mtj;
    m1b.gnd.mfCOM = m1b.c.mtj + m1b.mf.COM;
    m1b.gnd.mtpj = m1b.c.mtj + m1b.mf.mtpj;

    locSphere_3_r_new = locSphere_3_r - m1b.c.mtj';
    locSphere_4_r_new = locSphere_4_r - m1b.c.mtj';
    
    %% plot results in calcaneus ref frame
    % This is the ref frame attached to the calcn in both the old model and the
    % new model. The detailed model needs to be scaled.

    figure
    subplot(2,1,1)
    hold on
    grid on

    plot(m2.gnd.mtpj(1),m2.gnd.mtpj(2),'.b')
    plot(m2.gnd.mtj(1),m2.gnd.mtj(2),'.b')
    plot(m2.gnd.mfCOM(1),m2.gnd.mfCOM(2),'*b')
    plot(m2.gnd.cCOM(1),m2.gnd.cCOM(2),'*b')

    plot(m2.gnd_s.mtpj(1),m2.gnd_s.mtpj(2),'.k')
    plot(m2.gnd_s.mtj(1),m2.gnd_s.mtj(2),'.k')
    plot(m2.gnd_s.mfCOM(1),m2.gnd_s.mfCOM(2),'*k')
    plot(m2.gnd_s.cCOM(1),m2.gnd_s.cCOM(2),'*k')

    plot(m1b.gnd.mtpj(1),m1b.gnd.mtpj(2),'or')
    plot(m1b.c.mtj(1),m1b.c.mtj(2),'or')
    plot(m1b.gnd.mfCOM(1),m1b.gnd.mfCOM(2),'xr')
    plot(m1b.c.COM(1),m1b.c.COM(2),'xr')

    plot(m0.cmf.mtpj(1),m0.cmf.mtpj(2),'og')
    plot(m0.cmf.COM(1),m0.cmf.COM(2),'xg')
    plot(0,0,'og')
    title('side view (sagittal plane)')

    plot(locSphere_1_r(1),locSphere_1_r(2),'*c')
    plot(locSphere_2_r(1),locSphere_2_r(2),'*c')
    plot(locSphere_3_r(1),locSphere_3_r(2),'*c')
    plot(locSphere_4_r(1),locSphere_4_r(2),'*c')
    plot(m0.cmf.mtpj(1)+locSphere_5_r(1),m0.cmf.mtpj(2)+locSphere_5_r(2),'*c')
    plot(m0.cmf.mtpj(1)+locSphere_6_r(1),m0.cmf.mtpj(2)+locSphere_6_r(2),'*c')

    subplot(2,1,2)
    hold on
    grid on
    
    plot(m2.gnd.mtpj(1),m2.gnd.mtpj(3),'.b')
    plot(m2.gnd.mtj(1),m2.gnd.mtj(3),'.b')
    plot(m2.gnd.mfCOM(1),m2.gnd.mfCOM(3),'*b')
    plot(m2.gnd.cCOM(1),m2.gnd.cCOM(3),'*b')

    plot(m2.gnd_s.mtpj(1),m2.gnd_s.mtpj(3),'.k')
    plot(m2.gnd_s.mtj(1),m2.gnd_s.mtj(3),'.k')
    plot(m2.gnd_s.mfCOM(1),m2.gnd_s.mfCOM(3),'*k')
    plot(m2.gnd_s.cCOM(1),m2.gnd_s.cCOM(3),'*k')

    plot(m1b.gnd.mtpj(1),m1b.gnd.mtpj(3),'or')
    plot(m1b.c.mtj(1),m1b.c.mtj(3),'or')
    plot(m1b.gnd.mfCOM(1),m1b.gnd.mfCOM(3),'xr')
    plot(m1b.c.COM(1),m1b.c.COM(3),'xr')

    plot(m0.cmf.mtpj(1),m0.cmf.mtpj(3),'og')
    plot(m0.cmf.COM(1),m0.cmf.COM(3),'xg')
    plot(0,0,'og')
    title('bottom view')

    plot(locSphere_1_r(1),locSphere_1_r(3),'*c')
    plot(locSphere_2_r(1),locSphere_2_r(3),'*c')
    plot(locSphere_3_r(1),locSphere_3_r(3),'*c')
    plot(locSphere_4_r(1),locSphere_4_r(3),'*c')
    plot(m0.cmf.mtpj(1)+locSphere_5_r(1),m0.cmf.mtpj(3)+locSphere_5_r(3),'*c')
    plot(m0.cmf.mtpj(1)+locSphere_6_r(1),m0.cmf.mtpj(3)+locSphere_6_r(3),'*c')

    %% relating vectors to foot arch compression
    % based on DOI: 10.1038/srep19403

    a = m1b.c.mtj(1:2);
    b = m1b.mf.mtpj(1:2);

    l_0 = norm(a+b);
    h0 = -b(2); 
    % Difference in y-coordinate between calcn origin and m1a.f.mtpj 
    % results in 1° difference, so it is omitted.

    c0 = acos(h0/norm(a));
    d0 = acos(h0/norm(b));

    mt0 = (c0+d0)*180/pi;

    h1 = h0*0.8; 
    c1 = acos(h1/norm(a));
    d1 = acos(h1/norm(b));

    mt1 = (c1+d1)*180/pi;

    mt_bound = mt1 - mt0;
    % So 15° bound is sensible

    h2 = h0*0.87; 
    c2 = acos(h2/norm(a));
    d2 = acos(h2/norm(b));

    mt2 = (c2+d2)*180/pi;

    mt_2 = mt2 - mt0;
    
end


%%
% vars = x;
% vars(4) = vars(4)/2;
% vars(7) = vars(7)/2;
% vars(10) = vars(10)/2;
% 
% res = f_footmodel_add_mtj(vars,m0.cmf.m,m0.cmf.I,m0.cmf.COM,m2.c.m,m2.m.m,m2.f.m,m2.gnd.mtj2mfCOM,m2.gnd.mtj2cCOM,m1b.c.I,sf);
% res'