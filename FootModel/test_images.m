clear
close all
clc


scs = get(0,'ScreenSize');
fwide = [scs(3)/4*3, scs(4)/2];
fhigh = [scs(3)/2, scs(4)-120];
flong = [scs(3)/2, scs(4)/4];

fpos = [1,scs(4)/2+20;
        1,40;
        -scs(3)/2,40;
        -scs(3),40;];

pathmain = pwd;
[pathRepo,~,~]  = fileparts(pathmain);

folder = '\Figures';


file = 'PF_force_Caravaggi09.png';
pathRefImg = fullfile(pathRepo,folder,file);
img_F_PF = imread(pathRefImg);


x = [0; 100];
y = [0,1; 0,1];

figure('Position',[fpos(4,:),scs(3)/2, scs(4)/3]);
hold on
axis tight
hi1 = image([-13,112],flip([-0.16,1.75]),img_F_PF);
uistack(hi1,'bottom')
ylabel('Tension (%BW');
xlabel('Stance phase (%)');
plot(x,y) 
        