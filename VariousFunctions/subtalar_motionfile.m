[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);


N = 3;
nq.all = 33;

q_opt_GUI_GC = zeros(N,1+nq.all+2);




muscleNames = {'glut_med1_r','glut_med2_r','glut_med3_r',...
    'glut_min1_r','glut_min2_r','glut_min3_r','semimem_r',...
    'semiten_r','bifemlh_r','bifemsh_r','sar_r','add_long_r',...
    'add_brev_r','add_mag1_r','add_mag2_r','add_mag3_r','tfl_r',...
    'pect_r','grac_r','glut_max1_r','glut_max2_r','glut_max3_r',......
    'iliacus_r','psoas_r','quad_fem_r','gem_r','peri_r',...
    'rect_fem_r','vas_med_r','vas_int_r','vas_lat_r','med_gas_r',...
    'lat_gas_r','soleus_r','tib_post_r','flex_dig_r','flex_hal_r',...
    'tib_ant_r','per_brev_r','per_long_r','per_tert_r','ext_dig_r',...
    'ext_hal_r','ercspn_r','intobl_r','extobl_r','ercspn_l',...
    'intobl_l','extobl_l'};

NMuscle = length(muscleNames(1:end-3))*2;
Acts_GC = zeros(N,NMuscle);

pathOpenSim = [pathRepo,'/OpenSim'];
addpath(genpath(pathOpenSim));
JointAngle.labels = {'time','pelvis_tilt','pelvis_list',...
    'pelvis_rotation','pelvis_tx','pelvis_ty','pelvis_tz',...
    'hip_flexion_l','hip_adduction_l','hip_rotation_l',...
    'hip_flexion_r','hip_adduction_r','hip_rotation_r',...
    'knee_angle_l','knee_angle_r','ankle_angle_l','ankle_angle_r',...
    'subtalar_angle_l','subtalar_angle_r','mtj_angle_l','mtj_angle_r','mtp_angle_l','mtp_angle_r',...
    'lumbar_extension','lumbar_bending','lumbar_rotation',...
    'arm_flex_l','arm_add_l','arm_rot_l',...
    'arm_flex_r','arm_add_r','arm_rot_r',...
    'elbow_flex_l','elbow_flex_r',...
    'pro_sup_l','pro_sup_r'};

isubt = strcmp(JointAngle.labels,'subtalar_angle_r');
% Two gait cycles
% Joint angles
q_opt_GUI_GC_2 = [q_opt_GUI_GC;q_opt_GUI_GC];
q_opt_GUI_GC_2(:,1) = [0,1,2,3,4,5]';
q_opt_GUI_GC_2(1,isubt) = -20;
q_opt_GUI_GC_2(2,isubt) = -10;
q_opt_GUI_GC_2(3,isubt) = 0;
q_opt_GUI_GC_2(4,isubt) = 10;
q_opt_GUI_GC_2(5,isubt) = 20;
% Muscle activations (to have muscles turning red when activated).
Acts_GC_GUI = [Acts_GC;Acts_GC];
% Combine data joint angles and muscle activations
JointAngleMuscleAct.data = [q_opt_GUI_GC_2,Acts_GC_GUI];
% Get muscle labels
muscleNamesAll = cell(1,NMuscle);
for i = 1:NMuscle/2
    muscleNamesAll{i} = [muscleNames{i}(1:end-2),'_l'];
    muscleNamesAll{i+NMuscle/2} = [muscleNames{i}(1:end-2),'_r'];
end
% Combine labels joint angles and muscle activations
JointAngleMuscleAct.labels = JointAngle.labels;
for i = 1:NMuscle
    JointAngleMuscleAct.labels{i+size(q_opt_GUI_GC_2,2)} = ...
        [muscleNamesAll{i},'/activation'];
end
OutFolder = fullfile(pathRepo,'Results');
filenameJointAngles = fullfile(OutFolder,'subtalar.mot');
write_motionFile(JointAngleMuscleAct, filenameJointAngles);



