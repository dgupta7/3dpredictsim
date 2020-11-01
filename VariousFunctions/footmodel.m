function rest = footmodel(mc,mf,Ic,If,COMc,COMf,MTJ)
% variables:
% mc: mass calcaneus + midfoot
% mf: mass forefoot
% Ic: inertia (xx,yy,zz) calcaneus + midfoot
% If: inertia (xx,yy,zz) forefoot
% COMc: center of mass (x,y,z) calcaneus + midfoot
% COMf: center of mass (x,y,z) forefoot
% MTJ: position(x,y,z) of tarsometarsal joint in calcaneus reference frame
% 2 + 3*5 = 17dof

%% relations imposed by old, simpler footmodel
% parameters:
% m0: mass simplified calcaneus
% I0: inertia (xx,yy,zz) simplified calcaneus
% COM0: center of mass (x,y,z) simplified calcaneus
ankle0 = [0 -0.41085882914747662 0];  % tibia ref
subtalar0 = [-0.044572100000000003 -0.038339100000000001 -0.0072382799999999997];  %talus ref
mtp0 = [0.163409678774199 -0.00182784875586352 0.000987038328166303];   % calcn ref
tibia2mtp0 = ankle0 + subtalar0 + mtp0;


% combined center of mass remains the same
rest(1) = -COM0(1)*m0 + COMc(1)*mc + COMf(1)*mf;
rest(2) = -COM0(2)*m0 + COMc(2)*mc + COMf(2)*mf;
rest(3) = -COM0(3)*m0 + COMc(3)*mc + COMf(3)*mf;

% mass conservation
rest(4) = mc+mf-m0;

% combined inertia remains the same
rest(5) = Ic(1) + mc*sqrt(COMc(2)^2 + COMc(3)^2) + If(1) + mf*sqrt((COMf(2)+MTJ(2))^2 + (COMf(3)+MTJ(3))^2) - (I0(1) + m0*sqrt(COM0(2)^2 + COM0(3)^2));
rest(6) = Ic(2) + mc*sqrt(COMc(3)^2 + COMc(1)^2) + If(2) + mf*sqrt((COMf(3)+MTJ(3))^2 + (COMf(1)+MTJ(1))^2) - (I0(2) + m0*sqrt(COM0(3)^2 + COM0(1)^2));
rest(7) = Ic(3) + mc*sqrt(COMc(1)^2 + COMc(2)^2) + If(3) + mf*sqrt((COMf(1)+MTJ(1))^2 + (COMf(2)+MTJ(2))^2) - (I0(3) + m0*sqrt(COM0(1)^2 + COM0(2)^2));


%% relations based on detailed footmodel (more detailed than the model we want to obtain)
% parameters:
com_talus = [0.0139331 -0.00637782 -0.00622932];
com_calcn = [-0.0183999 -0.0127205 0.00168919];
com_midfoot = [0.0136264 -0.00541573 -0.00226789];
com_forefoot = [0.030362 -0.00692578 -0.000413788];

ankle = [0 -0.415 0];   % tibia ref
subtalar = [0.010347034099435 -0.031630326706915 0.00134647472860196];  %talus ref
chopart = [0.0221509217480428 0.00563051896570316 -0.000938103987425259];   % calcn ref
tarsometarsal = [0.0224650598180872 -0.0131925031448367 0.00457222612063733];   % midfoot ref
mtp = [0.0604862357021971 -0.0140878323268576 0.00286827055184947];     % frontfoot ref
m_calcn = 0.289923030444027;
m_midfoot = 0.139888239487345;
m_forefoot = 0.236187108315951;

% ratio of masses is the same
massratio = (m_calcn+m_midfoot)/m_forefoot;
rest(8) = mc/mf - massratio;

% center of mass of calcn and midfoot combined
tibia2mtp = ankle + subtalar + chopart + tarsometarsal + mtp; % tibia provides ref independent of foot model
tibia2com_calcn = ankle + subtalar + com_calcn;
tibia2com_midfoot = ankle + subtalar + chopart + com_midfoot;
tibia2com_cm = (tibia2com_calcn * m_calcn + tibia2com_midfoot * m_midfoot)*(m_calcn + m_midfoot);

rest(9) = (ankle0(1) + subtalar0(1) + COMc(1))/abs(tibia2mtp0) - tibia2com_cm(1)/abs(tibia2mtp);
rest(10) = (ankle0(2) + subtalar0(2) + COMc(2))/abs(tibia2mtp0) - tibia2com_cm(2)/abs(tibia2mtp);
rest(11) = (ankle0(3) + subtalar0(3) + COMc(3))/abs(tibia2mtp0) - tibia2com_cm(3)/abs(tibia2mtp);

% position tarsometarsal joint in tibia reference frame
tibia2tarsometarsal = ankle + subtalar + chopart + tarsometarsal;
rest(12) = (ankle0(1) + subtalar0(1) + MTJ(1))/abs(tibia2mtp0) - tibia2tarsometarsal(1)/abs(tibia2mtp);
rest(13) = (ankle0(2) + subtalar0(2) + MTJ(2))/abs(tibia2mtp0) - tibia2tarsometarsal(2)/abs(tibia2mtp);
rest(14) = (ankle0(3) + subtalar0(3) + MTJ(3))/abs(tibia2mtp0) - tibia2tarsometarsal(3)/abs(tibia2mtp);

%%
% 3 relations about inertias needed





end