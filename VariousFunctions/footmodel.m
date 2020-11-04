function F = footmodel(vars,m0,I0,COM0,mtp0,m_calcn,m_midfoot,m_forefoot,tarsometarsal2com_cm,com_forefoot,I_ff0,sf)
% variables:
mc  = vars(1);
mf  = vars(2);
Icx = vars(3);
Icy = vars(4);
Icz = vars(5);
Ifx = vars(6);
Ify = vars(7);
Ifz = vars(8);
COMcx = vars(9);
COMcy = vars(10);
COMcz = vars(11);
COMfx = vars(12);
COMfy = vars(13);
COMfz = vars(14);
MTJx = vars(15);
MTJy = vars(16);
MTJz = vars(17);

%% relations imposed by old, simpler footmodel
% parameters of new model
Ic = [Icx, Icy, Icz];
If = [Ifx, Ify, Ifz];
COMc = [COMcx, COMcy, COMcz];
COMf = [COMfx, COMfy, COMfz];
MTJ = [MTJx, MTJy, MTJz];
% derived parameters of new model
MTPJ = mtp0 - MTJ;
c2COMf = MTJ + COMf;
MTJ2COMc = -(MTJ - COMc);

% combined center of mass remains the same
F(1) = -COM0(1)*m0 + COMc(1)*mc + c2COMf(1)*mf;
F(2) = -COM0(2)*m0 + COMc(2)*mc + c2COMf(2)*mf;
F(3) = -COM0(3)*m0 + COMc(3)*mc + c2COMf(3)*mf;

% mass conservation
F(4) = mc+mf-m0;

% combined inertia remains the same
F(5) = Ic(1) + mc* (COMc(2)^2 + COMc(3)^2) + If(1) + mf* (c2COMf(2)^2 + c2COMf(3)^2) ... 
    - (I0(1) + m0* (COM0(2)^2 + COM0(3)^2));
F(6) = Ic(2) + mc* (COMc(3)^2 + COMc(1)^2) + If(2) + mf* (c2COMf(3)^2 + c2COMf(1)^2) ...
    - (I0(2) + m0* (COM0(3)^2 + COM0(1)^2));
F(7) = Ic(3) + mc* (COMc(1)^2 + COMc(2)^2) + If(3) + mf* (c2COMf(1)^2 + c2COMf(2)^2) ... 
    - (I0(3) + m0* (COM0(1)^2 + COM0(2)^2));


%% relations based on detailed footmodel (more detailed than the model we want to obtain)
% % parameters:
% com_calcn = [-0.0183999 -0.0127205 0.00168919];
% com_midfoot = [0.0136264 -0.00541573 -0.00226789];
% com_forefoot = [0.030362 -0.00692578 -0.000413788];
% chopart = [0.0221509217480428 0.00563051896570316 -0.000938103987425259];   % calcn ref
% tarsometarsal = [0.0224650598180872 -0.0131925031448367 0.00457222612063733];   % midfoot ref
% mtp = [0.0604862357021971 -0.0140878323268576 0.00286827055184947];     % frontfoot ref
% m_calcn = 0.289923030444027;
% m_midfoot = 0.139888239487345;
% m_forefoot = 0.236187108315951;


% ratio of masses is the same
massratio = (m_calcn+m_midfoot)/m_forefoot;
F(8) = mc/mf - massratio;

% center of mass of forefoot
F(9)  = COMf(1) - com_forefoot(1)*sf(1);
F(10) = COMf(2) - com_forefoot(2)*(sf(2));
F(11) = COMf(3) - com_forefoot(3)*sf(3);

% center of mass of calcn and midfoot combined 
% get params and vars in forefoot ref frame (tarsometarsal joint is origin)

% tarsometarsal2com_calcn = -(chopart + tarsometarsal - com_calcn);
% tarsometarsal2com_midfoot = -(tarsometarsal - com_midfoot);
% tarsometarsal2com_cm = (tarsometarsal2com_calcn * m_calcn + tarsometarsal2com_midfoot * m_midfoot)/(m_calcn + m_midfoot);

F(12) = MTJ2COMc(1) - tarsometarsal2com_cm(1)*sf(1);
F(13) = MTJ2COMc(2) - tarsometarsal2com_cm(2)*sf(2);
F(14) = MTJ2COMc(3) - tarsometarsal2com_cm(3)*sf(3);



%%
% % 3 relations about inertias needed: Open the forefoot of the detailed
% % footmodel in Fusion 360. Convert mesh to solid, and calculate its
% % inertias and mass for a unit mass density (assume constant mass distribution)
% 
% m_ff = 45.16;    % [g]
% I_ff = [2.351, 4.183, 2.718]*1e4;  % [g mm^2]
% 
% I_ff0 = I_ff/m_ff *1e-6;  % [kg m^2 /kg] mass invariant inertia

F(15) = If(1) - I_ff0(1) * mf;
F(16) = If(2) - I_ff0(2) * mf;
F(17) = If(3) - I_ff0(3) * mf;


end