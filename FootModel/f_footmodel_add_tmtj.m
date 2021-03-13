function F = footmodel(vars,m0,I0,COM0,m_calcn,m_midfoot,m_forefoot,tarsometatarsal2com_cm,com_forefoot,If,sf)
% This function discribes the relations between a foot model consisting or
% a calcaneus and toes, and a footmodel consisting of a calcaneus,
% metatarsi and toes.

% variables:
mc  = vars(1);
mf  = vars(2);
COMcx = vars(3);
COMcy = vars(4);
COMcz = vars(5);
COMfx = vars(6);
COMfy = vars(7);
COMfz = vars(8);
MTJx = vars(9);
MTJy = vars(10);
MTJz = vars(11);
Icx = vars(12);
Icy = vars(13);
Icz = vars(14);


%% relations imposed by old, simpler footmodel
% parameters of new model
Ic = [Icx, Icy, Icz];
COMc = [COMcx, COMcy, COMcz];
COMf = [COMfx, COMfy, COMfz];
MTJ = [MTJx, MTJy, MTJz];

% derived parameters of new model
c2COMf = MTJ + COMf;
MTJ2COMc = -(MTJ - COMc);

% combined center of mass remains the same
F(1) = -COM0(1)*m0 + COMc(1)*mc + c2COMf(1)*mf;
F(2) = -COM0(2)*m0 + COMc(2)*mc + c2COMf(2)*mf;
F(3) = -COM0(3)*m0 + COMc(3)*mc + c2COMf(3)*mf;

% mass conservation
F(4) = mc+mf-m0;

% combined inertia remains the same
F(12) = Ic(1) + mc* (COMc(2)^2 + COMc(3)^2) + If(1) + mf* (c2COMf(2)^2 + c2COMf(3)^2) ... 
    - (I0(1) + m0* (COM0(2)^2 + COM0(3)^2));
F(13) = Ic(2) + mc* (COMc(3)^2 + COMc(1)^2) + If(2) + mf* (c2COMf(3)^2 + c2COMf(1)^2) ...
    - (I0(2) + m0* (COM0(3)^2 + COM0(1)^2));
F(14) = Ic(3) + mc* (COMc(1)^2 + COMc(2)^2) + If(3) + mf* (c2COMf(1)^2 + c2COMf(2)^2) ... 
    - (I0(3) + m0* (COM0(1)^2 + COM0(2)^2));


%% relations based on detailed footmodel (more detailed than the model we want to obtain)
% ratio of masses is the same
massratio = (m_calcn+m_midfoot)/m_forefoot;
F(5) = mc/mf - massratio;

% center of mass of forefoot
F(6)  = COMf(1) - com_forefoot(1)*sf(1);
F(7) = COMf(2) - com_forefoot(2)*(sf(2));
F(8) = COMf(3) - com_forefoot(3)*sf(3);

% center of mass of calcn and midfoot combined 
% get params and vars in forefoot ref frame (tarsometatarsal joint is origin)

F(9) = MTJ2COMc(1) - tarsometatarsal2com_cm(1)*sf(1);
F(10) = MTJ2COMc(2) - tarsometatarsal2com_cm(2)*sf(2);
F(11) = MTJ2COMc(3) - tarsometatarsal2com_cm(3)*sf(3);


end