% Save the experimental data from subject 1 by Antoine Falisse in the same
% format as the data from subject 1 by K. Poggensee
clear
clc

[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
load([pathRepo '\ExperimentalData\ExperimentalData.mat']);

%%
Dat.Normal.gc.colheaders = ExperimentalData.Q.subject1.Qs.colheaders;

Dat.Normal.gc.Qall_mean = ExperimentalData.Q.subject1.Qs.mean.*(pi/180);
Dat.Normal.gc.Qall_std = ExperimentalData.Q.subject1.Qs.std.*(pi/180);

Dat.Normal.gc.Qdotall_mean = ExperimentalData.Q.subject1.Qdots.mean.*(pi/180);
Dat.Normal.gc.Qdotall_std = ExperimentalData.Q.subject1.Qdots.std.*(pi/180);

Dat.Normal.gc.Tall_bio_mean = ExperimentalData.Torques.subject1.mean;
Dat.Normal.gc.Tall_bio_std = ExperimentalData.Torques.subject1.std;
Dat.Normal.gc.Tall_mean = ExperimentalData.Torques.subject1.mean;
Dat.Normal.gc.Tall_std = ExperimentalData.Torques.subject1.std;

Dat.Normal.EMGheaders = ExperimentalData.EMG.subject1.colheaders;
Dat.Normal.gc.lowEMG_mean = ExperimentalData.EMG.subject1.mean;
Dat.Normal.gc.lowEMG_std = ExperimentalData.EMG.subject1.std;

Dat.Normal.gc.GRF.Fmean = ExperimentalData.GRFs.subject1.mean;
Dat.Normal.gc.GRF.Fstd = ExperimentalData.GRFs.subject1.std;

%%
P_mech = ExperimentalData.Q.subject1.Qdots.all.*(pi/180).*ExperimentalData.Torques.subject1.all;

P_mech_mean = mean(P_mech,3);
P_mech_std = std(P_mech,0,3);



%%
save([pathRepo '\Data\Fal_s1.mat','Dat']);