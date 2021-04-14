% This script derives polynomial models for the plantar fascia stiffness
% and the lumped midtarsal stiffness. The parameters are chosen such that
% they result in the observations reported by Ker et al[1] when they are
% used in the foot model as described here.

% [1] R. F. Ker, M. B. Bennett, S. R. Bibby, R. C. Kester and R. M. Alexander,
% „The spring in the arch of the human foot,” Nature, nr. 325, pp. 147-149, 1987. 
%
% Author: Lars D'Hondt (April 2021)

close all
clear
clc
pathmain        = pwd;
[pathRepo,~,~]  = fileparts(pathmain);

make_figures = 1;
sf = 0.90; % Ker's data is from 85kg male (shoe size 38 to 42 is 10% diff in foot length)

%% Fit load to torque relation
% get values from f_staticFootCompression_v4
load(fullfile(pathRepo,'Results','FootModel',...
    'Foot_3D_Fal_s1_mtj_subt1_v1_none_Ker1987_Q-30_30_F0_4000_WLv3_ls135.mat'),'R');

x = R.Fs_tib(R.Fs_tib<=3000)';
y = R.M_WL(2,(R.Fs_tib<=3000))';

[Torque,err1] = f_fitPolynomialFunction(x,y,3);

if make_figures
y1 = Torque(x);
figure
plot(x,y)
hold on
plot(x,y1,'--')
xlabel('vertical load (kN)')
ylabel('midtarsal torque (Nm)')
legend('actual','poly','Location','best')
title('Foot without plantar fascia')
end

%% Fit displacement to angle
x = R.l_fa(2,(R.Fs_tib<=3000))' - R.l_fa(2,1)';
y = R.Qs(2,(R.Fs_tib<=3000),R.jointfi.tmt.r)';

[angle,err2] = f_fitPolynomialFunction(x,y,5);

if make_figures
y1 = angle(x);
figure
plot(x,y)
hold on
plot(x,y1,'--')
xlabel('arch elongation (m)')
ylabel('midtarsal angle (rad)')
legend('actual','poly','Location','best')
title('Foot without plantar fascia')
end

%% reference figure
folder = '\Figures';
file = 'arch_stiffness_Ker87.png';
pathRefImg = fullfile(pathRepo,folder,file);
img = imread(pathRefImg);

h1 = figure('Position',[1000,200,600,660]);
hold on
h = image([1,9.65],flip([0,4]),img);
uistack(h,'bottom')
xlabel('horizontal elongation (mm)')
ylabel('vertical force (kN)')
title('Foot arch stiffness, as defined by Ker et al, 1987')

%% reference points 
% on graph c: plantar fascia and heel pad removed, everything else present

x = [1,    2,    3,    4,   5,   6,     7,    7.5,  8,    8.5, 8.8];
y = [0.05, 0.13, 0.27, 0.5, 0.77, 1.15, 1.65, 1.95, 2.35, 2.72, 3];

dl = linspace(0,9,20);
F_L = interp1(x,y,dl,'spline','extrap');

figure(h1)
plot(x,y,'o')
plot(dl,F_L)

x = x*sf; % elongation scales with size
y = y*sf^2; % force scales with surface

% transfer to rotational variables
q_mt_K = angle(x*1e-3);
T_mt_K = Torque(y*1e3);

h2 = figure;
plot(q_mt_K*180/pi,T_mt_K,'o','DisplayName','Ker 87')
hold on
xlabel('Angle (°)')
ylabel('Torque (Nm)')
title('Approximating arch stiffness (Ker, c) with midtarsal spring')
legend

qs = linspace(0,35,100)'*pi/180;
Ts = interp1(q_mt_K,T_mt_K,qs,'spline','extrap');
% plot(qs*180/pi,Ts,'DisplayName','Ker 87')

%% Fit curve in rotational variables
figure(h2)

q_mt = linspace(-40,40,500)'*pi/180;

% values for positive angles
xp = q_mt_K';
yp = T_mt_K';
% values for negative angles
xn = -flip(xp(1:end))*1; % stiffens half as fast
yn = -flip(yp(1:end)) + 2*Ts(1); % make them meet in the middle

x = [xn;xp];
y = [yn;yp];

plot(x*180/pi,y,'DisplayName','Both sides')
ylim([-50,50])

order = 11;
stop = 0;

while ~stop
    [T_pass_mtj,err3] = f_fitPolynomialFunction(x,y,order);
    M_li = T_pass_mtj(q_mt);
    plot(q_mt*180/pi,M_li,'DisplayName',[num2str(order) '^e order approx'])

    disp(err3)
    
    if err3 <= 0.1 || order >= 11
        stop = 1;
    else
        order = order+2;
    end
end


% get code string to paste into getPassiveMtjMomentWindlass_v3
[~,~,coeff] = f_fitPolynomialFunction(x,y,order);

% codeString = ['M_li = ' num2str(coeff(1),5)];
% for ii=1:length(coeff)-1
%     codeString = [codeString ' + ' num2str(coeff(ii+1),4) '*q_mt^' num2str(ii)];
% end
% codeString = [codeString ';'];
% disp(codeString)


% odd function with constant term
codeString = ['M_li = ' num2str(coeff(1),6)];
for ii=1:2:length(coeff)-1
    codeString = [codeString ' + ' num2str(coeff(ii+1),6) '*q_mt^' num2str(ii)];
end
codeString = [codeString ';'];
disp(codeString)

%% extract plantar fascia characteristic
load(fullfile(pathRepo,'Results','FootModel',...
    'Foot_3D_Fal_s1_mtj_subt1_v1_Ker1987_Ker1987_Q-30_30_F0_4000_WLv3_ls135.mat'),'R');


%% Fit load to torque relation
x = R.Fs_tib';
y = R.M_WL(2,:)';

[Torque1,err4] = f_fitPolynomialFunction(x,y,3);

if make_figures
y1 = Torque1(x);
figure
plot(x,y)
hold on
plot(x,y1,'--')
xlabel('vertical load (kN)')
ylabel('midtarsal torque (Nm)')
legend('actual','poly','Location','best')
title('Foot with plantar fascia')
end

%% Fit displacement to angle
x = R.l_fa(2,:)' - R.l_fa(2,1)';
y = R.Qs(2,:,R.jointfi.tmt.r)';

[angle1,err5] = f_fitPolynomialFunction(x,y,5);

if make_figures
y1 = angle1(x);
figure
plot(x,y)
hold on
plot(x,y1,'--')
xlabel('arch elongation (m)')
ylabel('midtarsal angle (rad)')
legend('actual','poly','Location','best')
title('Foot with plantar fascia')
end

%% reference points 
% on graph a: heel pad removed, everything else present

x = [0, 1,    2,    3,    4,    5, 6,   7,   7.5, 8, 8.5, 8.8];
y = [0, 0.05, 0.17, 0.37, 0.63, 1, 1.5, 2.1, 2.5, 3, 3.5, 3.8];

dl = linspace(0,9,20);
F_L = interp1(x,y,dl,'spline','extrap');

figure(h1)
plot(x,y,'o')
plot(dl,F_L)

% dl1 = [0.2,0.5];
% F_L1 = interp1(x,y,dl1,'spline');
% plot(dl1,F_L1,'d')
% 
% x = [0,dl1,x(2:end)];
% y = [0,F_L1,y(2:end)];

x = x*sf; % elongation scales with size
y = y*sf^2; % force scales with surface

% transfer to rotational variables
q_mt_K = angle1(x*1e-3)';
T_mt_K = Torque1(y*1e3)';

h3 = figure;
% plot(q_mt_K*180/pi,T_mt_K,'o','DisplayName','Total,Ker 87')
hold on
xlabel('Angle (°)')
ylabel('Torque (Nm)')
title('Approximating arch stiffness (Ker, a)')
legend

qs = linspace(0,35,100)'*pi/180;
Ts = interp1(q_mt_K,T_mt_K,qs,'spline','extrap');
plot(qs*180/pi,Ts,'DisplayName','Total, Ker 87')

M_li = T_pass_mtj(qs);
plot(qs*180/pi,M_li,'DisplayName','w/o PF, Ker 87')

M_PF = Ts - M_li;
plot(qs*180/pi,M_PF,'DisplayName','PF, Ker 87')

%% to plantar fascia variables
% same values as in getPassiveMtjMomentWindlass_v3 and f_staticFootCompression_v4
% windlass mechanism
calcnPF2mtj = 0.06695;
mtj2mttPF = 0.091714;
phi0 = 2.1274;
phi_K = phi0 + q_mt_K; % top angle of WL triangle
l_PF_fa_K = sqrt(calcnPF2mtj^2 + mtj2mttPF^2 - 2*calcnPF2mtj*mtj2mttPF*cos(phi_K)); % length of PF spanning arch
MA_PF_K = calcnPF2mtj*mtj2mttPF./l_PF_fa_K.*sin(phi_K); % moment arm of PF to mtj
q_mtp_K = -0.4751*q_mt_K; % mtp moves for toes to remain flat on ground
R_mtth = 7.5e-3; % average radius of the metatarsal head
l_PF_K = l_PF_fa_K + R_mtth*q_mtp_K;
lambda_PF_K = l_PF_K/R.S.PF_slack_length;
M_li_K = T_pass_mtj(q_mt_K);
M_PF_K = T_mt_K - M_li_K;
F_PF_K = -M_PF_K./MA_PF_K;


lambda_PF = linspace(1,lambda_PF_K(end)+0.01,100)';

%% Fit polynomial plantar fascia stiffness
h4 = figure;
hold on
plot(lambda_PF_K,F_PF_K,'o','DisplayName','PF, Ker 87')
legend('Location','northwest')
ylabel('Force (N)')
xlabel('Stretch ratio (-)')
title('Plantar fascia stiffness model')

x = [1; lambda_PF_K];
y = [0; F_PF_K];

% x = [lambda_PF_K];
% y = [F_PF_K];

order = 3;
stop = 0;
while ~stop
    [F_pass_PF,err5] = f_fitPolynomialFunction(x,y,order);
    F_PF = F_pass_PF(lambda_PF);
    plot(lambda_PF,F_PF,'DisplayName',[num2str(order) '^e order approx'])

    disp(err5)
    
    if err5 <= 0.1 || order >= 5
        stop = 1;
    else
        order = order+1;
    end
end

% get code string to paste into f_getPlantarFasciaStiffnessModelCasADiFunction.m
[~,~,coeff] = f_fitPolynomialFunction(x,y,order);
coeff(1) = coeff(1) - F_pass_PF(1); % l-ls = 0 -> F = 0

codeString = ['F_PF = ' num2str(coeff(1),10)];
for ii=1:length(coeff)-1
    codeString = [codeString ' + ' num2str(coeff(ii+1),10) '*lambda^' num2str(ii)];
end
codeString = [codeString ';'];
disp(codeString)















