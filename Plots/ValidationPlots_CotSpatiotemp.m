function [] = ValidationPlots_CotSpatiotemp(ResultsFiles,ReferenceFiles,DataFile)


load(DataFile,'Dat');

cot_diff_pas_dat = (Dat.Passive.COT - Dat.Normal.COT)/Dat.Normal.COT*100;
cot_diff_act_dat = (Dat.Active.COT - Dat.Normal.COT)/Dat.Normal.COT*100;

sf_diff_pas_dat = (Dat.Passive.Step.StrideFreq_mean - Dat.Normal.Step.StrideFreq_mean)/Dat.Normal.Step.StrideFreq_mean*100;
sf_diff_act_dat = (Dat.Active.Step.StrideFreq_mean - Dat.Normal.Step.StrideFreq_mean)/Dat.Normal.Step.StrideFreq_mean*100;

for i=1:3
    load(ResultsFiles{i},'R');
    if R.S.ExoBool
        if R.S.ExoScale
            cot_act = R.COT;
            sf_act = 1./(R.tf_step*2);
        else
            cot_pas = R.COT;
            sf_pas = 1./(R.tf_step*2);
        end
    else
        cot_nor = R.COT;
        sf_nor = 1./(R.tf_step*2);
    end
end

cot_diff_pas_res = (cot_pas - cot_nor)/cot_nor*100;
cot_diff_act_res = (cot_act - cot_nor)/cot_nor*100;

sf_diff_pas_res = (sf_pas - sf_nor)/sf_nor*100;
sf_diff_act_res = (sf_act - sf_nor)/sf_nor*100;

for i=1:3
    load(ReferenceFiles{i},'R');
    if R.S.ExoBool
        if R.S.ExoScale
            cot_act = R.COT;
            sf_act = 1./(R.tf_step*2);
        else
            cot_pas = R.COT;
            sf_pas = 1./(R.tf_step*2);
        end
    else
        cot_nor = R.COT;
        sf_nor = 1./(R.tf_step*2);
    end
end

cot_diff_pas_ref = (cot_pas - cot_nor)/cot_nor*100;
cot_diff_act_ref = (cot_act - cot_nor)/cot_nor*100;

sf_diff_pas_ref = (sf_pas - sf_nor)/sf_nor*100;
sf_diff_act_ref = (sf_act - sf_nor)/sf_nor*100;

figure(10)
subplot(1,2,1)
plot(0,cot_diff_pas_dat,'o')
hold on
grid on
plot(1,cot_diff_act_dat,'o')

plot(0,cot_diff_pas_ref,'*')
plot(1,cot_diff_act_ref,'*')

plot(0,cot_diff_pas_res,'.')
plot(1,cot_diff_act_res,'.')
ylabel('cot relative to normal')
xlim([-0.5,1.5])


subplot(1,2,2)
plot(0,sf_diff_pas_dat,'o')
hold on
grid on
plot(1,sf_diff_act_dat,'o')

plot(0,sf_diff_pas_ref,'*')
plot(1,sf_diff_act_ref,'*')

plot(0,sf_diff_pas_res,'.')
plot(1,sf_diff_act_res,'.')
ylabel('strid freq relative to normal')
xlim([-0.5,1.5])

end

