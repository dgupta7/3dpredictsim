% get muscle parameters
%------------------------


ModelPath = 'C:\Users\u0088756\Documents\FWO\Software\ExoSim\SimExo_3D\3dpredictsim\OpenSimModel';
OsimModel = fullfile(ModelPath,'Rajagopal2015.osim');
Mnames = DispMusclesOsimModel(osimModel);
muscleNames = Mnames(1:40);
MTparameters = ReadMuscleParameters(OsimModel,muscleNames);
save('MTparameters_Rajagopal2015.mat','MTparameters');
save('MuscleNames.mat','muscleNames');

