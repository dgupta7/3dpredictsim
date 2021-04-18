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

make_figures = 0;
sf = 0.9; % Ker's data is from 85kg male (shoe size 38 to 42 is 10% diff in foot length)
PF_slack_length = 0.144;


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

%% simulation reference

load(fullfile(pathRepo,'Results','FootModel',...
    'Foot_3D_Fal_s1_mtj_subt1_v1_none_Ker1987_Q-30_30_F0_4000_WLv3_ls137.mat'),'R');


j = find(R.Qs_mtp(:)==0);
js = find(R.failed(j,:)==0);
js = js(R.Fs_tib<=3000);

FL_sim = R.Fs_tib(js);
Q_mt_sim = squeeze(R.Qs(j,js,R.jointfi.tmt.r))*180/pi;
MA_sim_v = squeeze(R.metatarsi_or(j,js,[1,3]))-squeeze(R.talus_or(j,js,[1,3]));
L1_sim_v = squeeze(R.metatarsi_or(j,js,[1,3]))-squeeze(R.calcn_or(j,js,[1,3]));
L2_sim_v = -squeeze(R.metatarsi_or(j,js,[1,3]))+squeeze(R.toes_or(j,js,[1,3]));
for i=1:length(js)
    MA_sim(i) = norm(MA_sim_v(i,:));
    L1_sim(i) = norm(L1_sim_v(i,:));
    L2_sim(i) = norm(L2_sim_v(i,:));
end

[MA_p,err1] = f_fitPolynomialFunction(FL_sim',MA_sim',4);
MA_FL_p = MA_p(FL_sim');

if make_figures
    figure
    plot(FL_sim,MA_sim)
    hold on
    plot(FL_sim,MA_FL_p,'--')
    
    fc1 = figure;
    subplot(4,2,1)
    plot(Q_mt_sim,R.M_li(j,js))
    hold on
    % plot(q_mt_K*180/pi,T_mt_K,'o')
    xlabel('mtj angle (°)')
    ylabel('mtj torque (Nm)')

    subplot(4,2,2)
    plot(FL_sim,R.M_li(j,js))
    hold on
    % plot(FL,T_mt_K,'o')
    xlabel('load (N)')
    ylabel('mtj torque (Nm)')

    subplot(4,2,3)
    plot(Q_mt_sim,R.GRF_calcn(j,js,2))
    hold on
    % plot(q_mt_K*180/pi,R1,'o')
    xlabel('mtj angle (°)')
    ylabel('calcn GRF (N)')

    subplot(4,2,4)
    plot(FL_sim,R.GRF_calcn(j,js,2))
    hold on
    % plot(FL,R1,'o')
    xlabel('load (N)')
    ylabel('calcn GRF (N)')

    subplot(4,2,5)
    plot(Q_mt_sim,MA_sim)
    hold on
    % plot(q_mt_K*180/pi,MA_FL,'o')
    xlabel('mtj angle (°)')
    ylabel('lever arm (m)')

    subplot(4,2,6)
    plot(FL_sim,MA_sim)
    hold on
    % plot(FL,MA_FL,'o')
    xlabel('load (N)')
    ylabel('lever arm (m)')

    subplot(4,2,7)
    plot(FL_sim,L1_sim)
    hold on
    xlabel('load (N)')
    ylabel('l1 (m)')

    subplot(4,2,8)
    plot(FL_sim,L2_sim)
    hold on
    xlabel('load (N)')
    ylabel('l2 (m)')
    
end


%% Geometry constants
% run \FootModel\solveFootmodelParameters.m to get these values
% foot arch
calcn2mtj = 0.08207;
mtj2mtpj = 0.089638;
beta0 = 2.4935;
% windlass mechanism
calcnPF2mtj = 0.06695;
mtj2mttPF = 0.091714;
phi0 = 2.1274;

MA_FL0 = 0.0335; % lever arm of vertical external load around mtj in neutral position

%% reference points 
% on graph c: plantar fascia and heel pad removed, everything else present

x = [1,    2,    3,    4,   5,   6,     7,    7.5,  8,    8.5, 8.8];
y = [0.01, 0.13, 0.27, 0.5, 0.77, 1.15, 1.65, 1.95, 2.35, 2.72, 3];

% x1 = [x,9.5,9.8];
% y = interp1(x,y,x1,'spline','extrap');
% x = x1;

dl = linspace(0,9.8,100);
F_L = interp1(x,y,dl,'spline','extrap');

figure(h1)
plot(x,y,'o')
plot(dl,F_L)

dl_fa = x'*1e-3*sf;
FL = y'*1e3*sf^2;

%% From external load and displacement to mtj torque and rotation
l_fa_0 = sqrt(calcn2mtj^2 + mtj2mtpj^2 - 2*calcn2mtj*mtj2mtpj*cos(beta0)); % initial foot arch length
l_fa = l_fa_0 + dl_fa;
beta = acos( (l_fa.^2 - calcn2mtj^2 - mtj2mtpj^2)./(-2*calcn2mtj*mtj2mtpj) );
q_mt = beta - beta0; % mtj angle
h_fa_0 = calcn2mtj*mtj2mtpj./l_fa_0*sin(beta0); % initial foot arch height
ksi0 = asin(h_fa_0/calcn2mtj); % angle from horizontal to calcn2mtj
d = MA_FL0/cos(ksi0); % proj of MA_FL onto calcn2mtj
h_fa = calcn2mtj*mtj2mtpj./l_fa.*sin(beta);
ksi = asin(h_fa/calcn2mtj);
MA_FL = d*cos(ksi);
% MA_FL = MA_p(FL);
L1 = sqrt(calcn2mtj^2 - h_fa.^2);
L2 = l_fa - L1;
% L1 = L1 - 0.007;
R1 = FL.*(L2+MA_FL)./(L1+L2);
R2 = FL.*(L1-MA_FL)./(L1+L2);
TL = FL.*MA_FL - R1.*h_fa;


% MA_FL_p = MomentArm(q_mt);
% TL = FL.*MA_FL_p;


q_mt_K = q_mt;
T_mt_K = -TL;

% if make_figures
%     figure(fp1)
%     plot(q_mt,MA_FL,'DisplayName','planar geometry')
% end

%%
h2 = figure;
plot(q_mt_K*180/pi,T_mt_K,'o','DisplayName','Ker 87')
hold on
xlabel('Angle (°)')
ylabel('Torque (Nm)')
title('Approximating arch stiffness (Ker, c) with midtarsal spring')
legend
% ylim([-50,50])

q_max = q_mt_K(end);
qs = linspace(0,q_max,100)'*pi/180;
Ts = interp1(q_mt_K,T_mt_K,qs,'spline','extrap');
plot(qs*180/pi,Ts,'DisplayName','Ker 87')

%% Fit curve in rotational variables
corrf1 = 1;
corrf2 = 1;
corrf1 = interp1([q_mt_K(1),q_mt_K(end)],[1,0.55],q_mt_K,'spline','extrap');
corrf2 = interp1([q_mt_K(1),q_mt_K(4),q_mt_K(end)],[1.2,1.2,1],q_mt_K,'spline','extrap');

figure(h2)
q_mt = linspace(-1,1,500)'*q_max*1.5;
grid on

% values for positive angles
xp = q_mt_K;
yp = T_mt_K.*corrf1.*corrf2;
% values for negative angles
xn = -flip(xp(1:end))*1; % stiffens half as fast
yn = -flip(yp(1:end)) + 2*Ts(1); % make them meet in the middle

x = [xn;xp];
y = [yn;yp];

plot(x*180/pi,y,'DisplayName','Symmetric')


order = 9;
stop = 0;

while ~stop
    [T_pass_mtj,err2] = f_fitPolynomialFunction(x,y,order);
    M_li = T_pass_mtj(q_mt);
    plot(q_mt*180/pi,M_li,'DisplayName',[num2str(order) '^e order approx'])

    disp(err2)
    
    if err2 <= 0.1 || order >= 9
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
    codeString = [codeString ' + ' num2str(coeff(ii+1),8) '*q_mt^' num2str(ii)];
end
codeString = [codeString ';'];
disp(codeString)

%% compare with footSim results

if make_figures
    figure(fc1)
    subplot(4,2,1)
    plot(q_mt_K*180/pi,T_mt_K,'o')

    subplot(4,2,2)
    plot(FL,T_mt_K,'o')

    subplot(4,2,3)
    plot(q_mt_K*180/pi,R1,'o')

    subplot(4,2,4)
    plot(FL,R1,'o')

    subplot(4,2,5)
    plot(q_mt_K*180/pi,MA_FL,'o')

    subplot(4,2,6)
    plot(FL,MA_FL,'o')
    
    subplot(4,2,7)
    plot(FL,L1,'o')

    subplot(4,2,8)
    plot(FL,L2,'o')
    
    figure
    plot(FL_sim,Q_mt_sim)
    hold on
    plot(FL,q_mt_K*180/pi,'o')
    
end

%% reference points 
% on graph a: heel pad removed, everything else present

x = [0, 1,    2,    3,    4,    5, 6,   7,   7.5, 8, 8.5, 8.8];
y = [0, 0.05, 0.17, 0.37, 0.63, 1, 1.5, 2.1, 2.5, 3, 3.5, 3.8];

% y1 = [y,5];
% x = interp1(y,x,y1,'spline','extrap');
% y = y1;

dl = linspace(0,x(end),100);
F_L = interp1(x,y,dl,'spline','extrap');

figure(h1)
plot(x,y,'o')
plot(dl,F_L)

dl_fa = x'*1e-3*sf;
FL = y'*1e3*sf^2;

%% From external load and displacement to mtj torque and rotation
l_fa_0 = sqrt(calcn2mtj^2 + mtj2mtpj^2 - 2*calcn2mtj*mtj2mtpj*cos(beta0)); % initial foot arch length
l_fa = l_fa_0 + dl_fa;
beta = acos( (l_fa.^2 - calcn2mtj^2 - mtj2mtpj^2)./(-2*calcn2mtj*mtj2mtpj) );
q_mt = beta - beta0; % mtj angle
h_fa_0 = calcn2mtj*mtj2mtpj./l_fa_0*sin(beta0); % initial foot arch height
ksi0 = asin(h_fa_0/calcn2mtj); % angle from horizontal to calcn2mtj
d = MA_FL0/cos(ksi0); % proj of MA_FL onto calcn2mtj
h_fa = calcn2mtj*mtj2mtpj./l_fa.*sin(beta);
ksi = asin(h_fa/calcn2mtj);
MA_FL = d*cos(ksi);
% MA_FL = MA_p(FL);
L1 = sqrt(calcn2mtj^2 - h_fa.^2);
L2 = l_fa - L1;
% L1 = L1 - 0.007;
R1 = FL.*(L2+MA_FL)./(L1+L2);
R2 = FL.*(L1-MA_FL)./(L1+L2);
TL = FL.*MA_FL - R1.*h_fa;

q_mt_K = q_mt;
T_mt_K = -TL;


%%
h3 = figure;
% plot(q_mt_K*180/pi,T_mt_K,'o','DisplayName','Total,Ker 87')
hold on
xlabel('Angle (°)')
ylabel('Torque (Nm)')
title('Approximating arch stiffness (Ker, a)')
legend
% ylim([-50,50])

qs = linspace(0,q_mt_K(end),100)';
Ts = interp1(q_mt_K,T_mt_K,qs,'spline','extrap');
plot(qs*180/pi,Ts,'DisplayName','Total, Ker 87')

M_li = T_pass_mtj(qs);
plot(qs*180/pi,M_li,'DisplayName','w/o PF, Ker 87')

M_PF = Ts - M_li;
plot(qs*180/pi,M_PF,'DisplayName','PF, Ker 87')
grid on

%% to plantar fascia variables
phi_K = phi0 + q_mt_K; % top angle of WL triangle
l_PF_fa_K = sqrt(calcnPF2mtj^2 + mtj2mttPF^2 - 2*calcnPF2mtj*mtj2mttPF*cos(phi_K)); 
        % length of PF spanning arch
MA_PF_K = calcnPF2mtj*mtj2mttPF./l_PF_fa_K.*sin(phi_K); % moment arm of PF to mtj
q_mtp_K = -0.4751*q_mt_K; % mtp moves for toes to remain flat on ground
R_mtth = 7.5e-3; % average radius of the metatarsal head
l_PF_K = l_PF_fa_K + R_mtth*(pi/2+q_mtp_K)-0.0017;
lambda_PF_K = l_PF_K/PF_slack_length;
M_li_K = T_pass_mtj(q_mt_K)-1;
M_PF_K = T_mt_K - M_li_K;
F_PF_K = -M_PF_K./MA_PF_K;

% some adjusting
F_PF_K = F_PF_K + 50;
F_PF_K(1) = F_PF_K(1) -75;
F_PF_K(2) = F_PF_K(2)/3;

% figure
% plot(q_mt_K*180/pi,lambda_PF_K)
% ylabel('Stretch ratio (-)')
% xlabel('Angle (°)')

lambda_PF = linspace(1,1.2,1000)';

%% Fit polynomial plantar fascia stiffness

corrf3 = 1;
corrf3 = interp1([q_mt_K(1),q_mt_K(4),q_mt_K(7),q_mt_K(9),q_mt_K(end)],[0.9,1,1,0.7,0.5],q_mt_K,'linear','extrap');

h4 = figure;
hold on
plot(lambda_PF_K,F_PF_K,'o','DisplayName','PF, Ker 87')
plot(lambda_PF_K,F_PF_K.*corrf3,'o','DisplayName','PF, Ker 87, adjusted')
legend('Location','northwest')
ylabel('Force (N)')
xlabel('Stretch ratio (-)')
title('Plantar fascia stiffness model')
grid on
ylim([-100,2e3]);
F_PF_K = F_PF_K.*corrf3;
x = [1;1.005;1.01; lambda_PF_K];
y = [0;1;3; F_PF_K];

% x = [lambda_PF_K];
% y = [F_PF_K];


order = 6;
stop = 0;
while ~stop
    [F_pass_PF,err3] = f_fitPolynomialFunction(x,y,order);
    F_PF = F_pass_PF(lambda_PF);
    plot(lambda_PF,F_PF,'DisplayName',[num2str(order) '^e order approx'])

    disp(err3)
    
    if err3 <= 0.1 || order >= 6
        stop = 1;
    else
        order = order+1;
    end
end

% get code string to paste into f_getPlantarFasciaStiffnessModelCasADiFunction.m
[~,~,coeff] = f_fitPolynomialFunction(x,y,order);
dFdl = zeros(length(lambda_PF),1);
for i=2:length(coeff)
    dFdl(:,1) = dFdl(:,1) + coeff(i).*lambda_PF.^(i-2)*(i-1);
end

figure
plot(lambda_PF,dFdl)
grid on
xlabel('\lambda = stretch ratio (-)')
ylabel('dF/d\lambda')
title('must be positive, so stiffness is strictly increasing')
ylim([-100,2e3]);

coeff(1) = coeff(1) - F_pass_PF(1); % l-ls = 0 -> F = 0

codeString = ['F_PF = ' num2str(coeff(1),11)];
for ii=1:length(coeff)-1
    codeString = [codeString ' + ' num2str(coeff(ii+1),12) '*lambda^' num2str(ii)];
end
codeString = [codeString ';'];
disp(codeString)



%%

%% extract plantar fascia characteristic
% load(fullfile(pathRepo,'Results','FootModel',...
%     'Foot_3D_Fal_s1_mtj_subt1_v1_Ker1987_Ker1987_Q-30_30_F0_4000_WLv3_ls137.mat'),'R');

% load(fullfile(pathRepo,'Results','FootModel',...
%     'Foot_3D_Fal_s1_mtj_subt1_v1_none_Ker1987_Q-30_30_F0_4000_WLv3_ls137.mat'),'R');
% 
% 
% j = find(R.Qs_mtp(:)==0);
% js = find(R.failed(j,:)==0);
% 
% FL_sim = R.Fs_tib(js);
% Q_mt_sim = squeeze(R.Qs(j,js,R.jointfi.tmt.r))*180/pi;
% 
% figure
% subplot(3,2,1)
% plot(Q_mt_sim,R.M_WL(j,js))
% hold on
% plot(q_mt_K*180/pi,T_mt_K,'o')
% 
% subplot(3,2,2)
% plot(FL_sim,R.M_WL(j,js))
% hold on
% plot(FL,T_mt_K,'o')

% subplot(3,2,3)
% plot(Q_mt_sim,R.M_PF(j,js))
% hold on
% plot(q_mt_K*180/pi,M_PF_K,'o')
% 
% subplot(3,2,4)
% plot(FL_sim,R.M_PF(j,js))
% hold on
% plot(FL,M_PF_K,'o')

% subplot(3,2,5)
% plot(Q_mt_sim,R.M_li(j,js))
% hold on
% plot(q_mt_K*180/pi,M_li_K,'o')
% 
% subplot(3,2,6)
% plot(FL_sim,R.M_li(j,js))
% hold on
% plot(FL,M_li_K,'o')


% figure
% plot(FL_sim,R.GRF_calcn(j,js,2))
% hold on
% plot(FL,R1,'o')




