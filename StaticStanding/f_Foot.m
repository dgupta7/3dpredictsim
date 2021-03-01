function [res] = f_Foot(vars,F,PassT_Foot,jointi,Qmtp,Ftibia)

nqall      = 10;

%%
Qs = zeros(nqall,1);
Qdots = zeros(nqall,1);
A = zeros(nqall,1);

Qs(1:9) = vars(:);
Qs(10) = Qmtp;

QsQdots_nsc = zeros(nqall*2,1);
QsQdots_nsc(1:2:end,:) = Qs;
QsQdots_nsc(2:2:end,:) = Qdots;

%%
Qs_pass = Qs([(jointi.ankle.r),(jointi.subt.r),(jointi.tmt.r),(jointi.mtp.r)]);
Tau_pass = PassT_Foot(Qs_pass);

%%
[Tj] = full(F([QsQdots_nsc(:,1);A(:,1);0]));

GRF_calcn(:) = Tj(jointi.calcn_GRF);
GRF_metatarsi(:) = Tj(jointi.metatarsi_GRF);
m_foot = 2.78+0.077+0.625+0.344+0.168;

toes_or = Tj(jointi.toes_or);
calcn_or = Tj(jointi.calcn_or);

res(1) = Qs(1)^2 + Qs(2)^2; % tibia orientation
res(2) = Qs(4)^2 + Qs(6)^2; % tibia position x,z
% res(5) = (tanh((Qs(5)-0.3)*100)-1)*100; % tibia position y > 0.3
% res(3) = (Ftibia + m_foot*9.81) - (GRF_calcn(2)+GRF_metatarsi(2)); % force equilibrium y
res(3) = Tj(jointi.tibia.ty,1)+Ftibia;
% res(4) = Tj(1)^2 + Tj(2)^2 + Tj(3)^2;
res(4) = Qs(3)*180/pi + 10;
res(5) = toes_or(2) - (0.02+0.01);
res(6) = calcn_or(2) - (0.017+0.01);
res(7) = Tj(jointi.ankle.r,1)-Tau_pass(1); % ankle
res(8) = Tj(jointi.subt.r,1)-Tau_pass(2); % subt
res(9) = Tj(jointi.tmt.r,1) - Tau_pass(3); % tmt
% mtp is fixed in position by external forces



end