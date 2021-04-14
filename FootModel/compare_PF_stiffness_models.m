clear
clc
AddCasadiPaths();

ls = 0.137;
l = linspace(ls,ls+0.03,1000);

% PF_stiffness = {'linear','tanh','Gefen2001','Cheng2008','Barrett2018','Natali2010'};
PF_stiffness = {'Cheng2008','Gefen2001','Ker1987','Natali2010','linear'};

CsV1 = hsv(numel(PF_stiffness));
A0 = 49.7;
dl = (l-ls)*1000;
figure
hold on
for i=1:numel(PF_stiffness)
    f_PF_stiffness = f_getPlantarFasciaStiffnessModelCasADiFunction(PF_stiffness{i},'ls',ls);
    for j=1:1000
        F_PF(i,j) = full(f_PF_stiffness(l(j)));
    end
    plot((l-ls)*1000,F_PF(i,:),'DisplayName',PF_stiffness{i},'color',CsV1(i,:))
end

legend('Location','best')
xlabel('Elongation (mm)')
ylabel('Force (N)')
title('Plantar fascia stiffness models')
ylim([0,2000])
% xlim([0,20])

%%

% sf = 3;
% Fsf = F_PF(3,:).*sf;
% plot((l-ls)*1000,Fsf,'DisplayName',['Gefen x' num2str(sf)])

