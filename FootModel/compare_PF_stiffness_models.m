
ls = 0.16;
l = linspace(ls,ls+0.03,1000);

% PF_stiffness = {'linear','tanh','sqr','exp','Gefen2001','Cheng2008','Barrett2018','Natali2010'};
PF_stiffness = {'Cheng2008','Barrett2018','Gefen2001','Natali2010','linear'};

CsV1 = hsv(numel(PF_stiffness));
A0 = 58.6;
dl = (l-ls)*1000;
figure
hold on
for i=1:numel(PF_stiffness)
    f_PF_stiffness = f_getPlantarFasciaStiffnessModelCasADiFunction(PF_stiffness{i});
    for j=1:1000
        [F,L] = f_PF_stiffness(l(j));
        F_PF(i,j) = full(F);
    end
    plot((l/ls-1)*100,F_PF(i,:)/A0,'DisplayName',PF_stiffness{i},'color',CsV1(i,:))
end

legend('Location','best')
% xlabel('\Deltal (mm)')
xlabel('Nominal strain (%)')
ylabel('Nominal stress (MPa)')
title('Plantar fascia stiffness models')
ylim([0,60])
% xlim([0,10])




