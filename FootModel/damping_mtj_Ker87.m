% This script derives a viscous friction coefficient for the midtarsal joint, 
% based on observations reported by Ker et al[1].

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

h1 = figure('Position',[1000,100,800,880]);
hold on
h = image([1,9.65],flip([0,4]),img);
uistack(h,'bottom')
xlabel('horizontal elongation (mm)')
ylabel('vertical force (kN)')
title('Foot arch stiffness, as defined by Ker et al, 1987')

f = 0.2; % Hz (see article)

%% reference points 
% on graph c: plantar fascia and heel pad removed, everything else present
% increasing
x1 = [1,    2,    3,    4,   5,    6,    7,    7.5,  8,    8.5,  8.8,  8.9];
y1 = [0.01, 0.13, 0.27, 0.5, 0.77, 1.15, 1.65, 1.95, 2.35, 2.72, 2.98, 3];

d1 = linspace(1,8.9,1000);
F1 = interp1(x1,y1,d1,'spline','extrap');

figure(h1)
plot(x1,y1,'-o')

% decreasing
x2 = [1,    2,   3,    4,    5,   6,    7,    7.5,  8,    8.5,  8.8, 8.9];
y2 = [0.01, 0.1, 0.19, 0.35, 0.6, 0.91, 1.35, 1.63, 1.98, 2.45, 2.8, 3];

d2 = linspace(1,8.9,1000);
F2 = interp1(x2,y2,d2,'spline','extrap');

figure(h1)
plot(x2,y2,'-o')


%% reference points 
% on graph f
% increasing
x3 = [2,    3,    4,    5,   6,    7,    8,    8.5,  8.9];
y3 = [0.01, 0.05, 0.07, 0.1, 0.17, 0.25, 0.36, 0.43, 0.48];

d3 = linspace(2,8.9,1000);
F3 = interp1(x3,y3,d3,'spline','extrap');

figure(h1)
plot(x3,y3,'-o')

% increasing
x4 = [2,    3,    4,    5,    6,    7,    8,    8.5,  8.9];
y4 = [0.01, 0.04, 0.05, 0.08, 0.13, 0.20, 0.30, 0.38, 0.48];

d4 = linspace(2,8.9,1000);
F4 = interp1(x4,y4,d4,'spline','extrap');

figure(h1)
plot(x4,y4,'-o')


%% energy dissipated
% graph c
W_in_c = trapz(d1*1e-3*0.9,F1*1e3*0.9^2);
W_out_c = trapz(d2*1e-3*0.9,F2*1e3*0.9^2);
E_c = W_in_c - W_out_c;

% graph f
W_in_f = trapz(d3*1e-3*0.9,F3*1e3*0.9^2);
W_out_f = trapz(d4*1e-3*0.9,F4*1e3*0.9^2);
E_f = W_in_f - W_out_f;

% by midtarsal joint
E = E_c - E_f;

%%

Fmean = (F1+F2)/2*1e3*0.9^2; 

dL = d1*1e-3*0.9;
N = 1000;
t = linspace(0,1/f,N);

F = Fmean(end)*(1 - cos(2*pi*f*t))/2;

ts1 = interp1(F(1:N/2),t(1:N/2),Fmean,'spline','extrap');
ts2 = interp1(F(N/2+1:end),t(N/2+1:end),flip(Fmean(1:end-1)),'spline','extrap');

ts = [0,ts1,ts2,t(end)];

% figure
% plot(t,F)

%%
load(fullfile(pathRepo,'Results','FootModel',...
    'Foot_3D_Fal_s1_mtj_subt1_v5_none_Ker1987_Q-30_30_F0_2500_WLv3_ls150_sb1.mat'),'R');


j = find(R.Qs_mtp(:)==0);
js = find(R.failed(j,:)==0);
js = js(R.Fs_tib<=3000);

dl_fa_sim = R.l_fa(j,js)-R.L0;
FL_sim = R.Fs_tib(js);
Q_mt_sim = squeeze(R.Qs(j,js,R.jointfi.tmt.r));

qs = interp1(dl_fa_sim,Q_mt_sim,dL,'spline','extrap');


%%

Qs = [qs(1),qs,flip(qs(1:end-1)),qs(1)];
Qdots = zeros(size(Qs));

for i=2:length(Qs)
   Qdots(i) = (Qs(i)-Qs(i-1)) / (ts(i)-ts(i-1));
end

Qdots(2) = interp1(ts([1,3:10]),Qdots([1,3:10]),ts(2),'spline');

figure
subplot(211)
plot(ts,Qs)
hold on
plot(t,F*max(Qs)/max(F),'--')

subplot(212)
plot(ts,Qdots)

%%

dMT = 1;

P = dMT*Qdots.^2; % dissipated power
E_d = trapz(ts,P);

dMT = E/E_d;













