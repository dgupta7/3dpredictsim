function F = f_footmodel_add_mtj(vars,m0,I0,COM0,m_calcn,m_midfoot,m_forefoot,mtj2com_midforefoot,mtj2com_calcn,Ic,sf)
% This function discribes the relations between a foot model consisting or
% a calcaneus and toes, and a footmodel consisting of a calcaneus,
% midfoot and toes.

% variables:
mc  = vars(1);
mmf  = vars(2);
COMcx = vars(3);
COMcy = vars(4);
COMcz = vars(5);
COMmfx = vars(6);
COMmfy = vars(7);
COMmfz = vars(8);
MTJx = vars(9);
MTJy = vars(10);
MTJz = vars(11);
Imfx = vars(12);
Imfy = vars(13);
Imfz = vars(14);


%% relations imposed by old, simpler footmodel
% parameters of new model
Imf = [Imfx, Imfy, Imfz];
COMc = [COMcx, COMcy, COMcz];
COMmf = [COMmfx, COMmfy, COMmfz];
MTJ = [MTJx, MTJy, MTJz];

% derived parameters of new model
c2COMmf = MTJ + COMmf;
MTJ2COMc = -MTJ + COMc;

% combined center of mass remains the same
F(1) = -COM0(1)*m0 + COMc(1)*mc + c2COMmf(1)*mmf;
F(2) = -COM0(2)*m0 + COMc(2)*mc + c2COMmf(2)*mmf;
F(3) = -COM0(3)*m0 + COMc(3)*mc + c2COMmf(3)*mmf;

% mass conservation
F(4) = mc+mmf-m0;

% combined inertia remains the same
F(12) = Ic(1) + mc* (COMc(2)^2 + COMc(3)^2) + Imf(1) + mmf* (c2COMmf(2)^2 + c2COMmf(3)^2) ... 
    - (I0(1) + m0* (COM0(2)^2 + COM0(3)^2));
F(13) = Ic(2) + mc* (COMc(3)^2 + COMc(1)^2) + Imf(2) + mmf* (c2COMmf(3)^2 + c2COMmf(1)^2) ...
    - (I0(2) + m0* (COM0(3)^2 + COM0(1)^2));
F(14) = Ic(3) + mc* (COMc(1)^2 + COMc(2)^2) + Imf(3) + mmf* (c2COMmf(1)^2 + c2COMmf(2)^2) ... 
    - (I0(3) + m0* (COM0(1)^2 + COM0(2)^2));


%% relations based on detailed footmodel (more detailed than the model we want to obtain)
% ratio of masses is the same
massratio = (m_calcn)/(m_forefoot+m_midfoot);
F(5) = mc/mmf - massratio;

% center of mass of forefoot and midfoot combined
% get params and vars in midfoot ref frame (midtarsal joint is origin)
F(6) = COMmf(1) - mtj2com_midforefoot(1)*sf(1);
F(7) = COMmf(2) - mtj2com_midforefoot(2)*sf(2);
F(8) = COMmf(3) - mtj2com_midforefoot(3)*sf(3);

% center of mass of calcn
% get params and vars in midfoot ref frame (midtarsal joint is origin)
F(9)  = MTJ2COMc(1) - mtj2com_calcn(1)*sf(1);
F(10) = MTJ2COMc(2) - mtj2com_calcn(2)*sf(2);
F(11) = MTJ2COMc(3) - mtj2com_calcn(3)*sf(3);


end