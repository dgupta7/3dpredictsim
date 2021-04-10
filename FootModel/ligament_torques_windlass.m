% Calculate midtarsal stiffness from Welte et al and Natali et al

clear
clc
AddCasadiPaths();


%% foot arch
calcn2mtj = 0.08207;
mtj2mtpj = 0.089638;
beta0 = 2.4935;

beta = beta0; % top angle of foot arch triangle
l_fa_0 = sqrt(calcn2mtj^2 + mtj2mtpj^2 - 2*calcn2mtj*mtj2mtpj*cos(beta)); % foot arch length

l_fa_max = l_fa_0/0.98;
l_fa_min = l_fa_max*0.96;

beta_max = acos( (l_fa_max^2 - calcn2mtj^2 - mtj2mtpj^2)/(-2*calcn2mtj*mtj2mtpj) );
q_mt_max = beta_max - beta0;

beta_min = acos( (l_fa_min^2 - calcn2mtj^2 - mtj2mtpj^2)/(-2*calcn2mtj*mtj2mtpj) );
q_mt_min = beta_min - beta0;

%
q_mtp_th = [-30,0,30]*pi/180;
q_mt_th = [q_mt_max,0,q_mt_min];
l_th = [l_fa_max,l_fa_0,l_fa_min];

l_1 = interp1(q_mtp_th,l_th,[-20,-10,0,10,20]*pi/180);
q_mtp = [-30:10:30]'*pi/180;
l = [l_fa_max,l_1,l_fa_min]';
q_mt = acos( (l.^2 - calcn2mtj^2 - mtj2mtpj^2)/(-2*calcn2mtj*mtj2mtpj) ) - beta0;



%% plantar fascia model
ls = 0.139;
PF_stiffness = {'Natali2010'};
f_PF_stiffness = f_getPlantarFasciaStiffnessModelCasADiFunction(PF_stiffness,'ls',ls);

S.sf_PF = 1;
S.sf_li = 0;
S.dMT = 0;
S.MT_li_nonl = 1;
S.kMT_li = 90;


for j=1:length(q_mtp)
    [M0,~,M1,M2,~,~,~,~,~,h1,l1] = ...
        getPassiveMtjMomentWindlass_v3(q_mt(j),0,q_mtp(j),f_PF_stiffness,S);
    M_PF(j) = M1;
    M_li(j) = M2;
    M(j) = M0;
    l_fa(j) = l1;
    h_fa(j) = h1;
 
end

f1 = figure;
plot(q_mt*180/pi,-M_PF,'o')
hold on
grid on
ylim([-60,60])

%%
qs_mt = linspace(-20,20,1000)'*pi/180;
k = 500;
plot(qs_mt*180/pi,-k*qs_mt)

%%
M_li = -2*exp(10*(qs_mt-5*pi/180)) + 2*exp(-15*(qs_mt+5*pi/180));

figure(f1)
plot(qs_mt*180/pi,M_li)



