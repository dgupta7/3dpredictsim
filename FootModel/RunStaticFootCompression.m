
clear
clc

%%
S.subject = 'Fal_s1';
S.Foot.Model = 'mtj';
S.Foot.Scaling = 'custom';

S.Foot.PF_stiffness = 'Gefen2002';
S.Foot.PF_sf = 1; 
S.Foot.PF_slack_length = 0.147;

S.Foot.MT_li_nonl = 1;
% S.Foot.mtj_stiffness = 'MG_exp';
S.Foot.mtj_stiffness = 'MG_poly';
% S.Foot.mtj_stiffness = 'Gefen2002';

%%
Qs_mtp = [0]*pi/180;
Fs_tib = [0];

% mtp angles to be considered
% Qs_mtp = [-45:15:45]*pi/180;
% Qs_mtp = [-30:30:30]*pi/180;
% Qs_mtp = [-20:5:30]*pi/180;
% Qs_mtp = [0:5:30]*pi/180;

% vertical forces on knee
% Fs_tib = [0:100:1000];
% Fs_tib = [0,100,320,640,960];
% Fs_tib = [0:0.1:1.5]*10*round(BW/10);
Fs_tib = [0:50:300,400:100:1000,1200:200:3000];
% Fs_tib = [0:50:300,400:100:1000,1250:250:2700];
% Fs_tib = [0:50:1000,1100:100:6000];



subtR = 1; % reduce subtalar mobility




%%


S.Foot.PF_stiffness = 'none';
R = f_staticFootCompression_v6(S,Qs_mtp,Fs_tib,subtR);
Results{1} = R;

S.Foot.PF_stiffness = 'Gefen2002';
R = f_staticFootCompression_v6(S,Qs_mtp,Fs_tib,subtR);
Results{2} = R;


%%
% call plot function
nrf = length(Results);
CsV = hsv(nrf);
for i=1:nrf
    R = Results{i};
    if i==1
        h = PlotResults_FootSim(R,CsV(i,:),0,0);
    else
        PlotResults_FootSim(R,CsV(i,:),h,0);
    end
end














