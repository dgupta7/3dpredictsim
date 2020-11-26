% This function makes default IK files that are compatible with
% PredSim_gait92_tmt

clear all
close all
clc

[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
cd([pathRepo '/OpenSimModel']);

%
load('IK_Bounds_Default.mat','Qsall');

Qsall.data(:,22:34) = Qsall.data(:,20:32);

for i=34:-1:22
    for j=1:6
    Qsall.textdata{j,i} = Qsall.textdata{j,i-2};
    end
    Qsall.colheaders{1,i} = Qsall.colheaders{1,i-2};    
end

Qsall.textdata{6,20} = 'tmt_angle_l';
Qsall.textdata{6,21} = 'tmt_angle_r';

Qsall.colheaders{1,20} = 'tmt_angle_l';
Qsall.colheaders{1,21} = 'tmt_angle_r';

save('IK_Bounds_Default_tmt.mat','Qsall');

%
load('IK_Guess_Default.mat','Qsall');

Qsall.data(:,22:34) = Qsall.data(:,20:32);

for i=34:-1:22
    for j=1:6
    Qsall.textdata{j,i} = Qsall.textdata{j,i-2};
    end
    Qsall.colheaders{1,i} = Qsall.colheaders{1,i-2};    
end

Qsall.textdata{6,20} = 'tmt_angle_l';
Qsall.textdata{6,21} = 'tmt_angle_r';

Qsall.colheaders{1,20} = 'tmt_angle_l';
Qsall.colheaders{1,21} = 'tmt_angle_r';

save('IK_Guess_Default_tmt.mat','Qsall');


