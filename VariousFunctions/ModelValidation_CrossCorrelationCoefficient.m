function [ccc] = ModelValidation_CrossCorrelationCoefficient(pathSimResults, pathExpData)
% returns a struct with information about the cross-correlation coefficient
% of the simulation results and the mean measurement data


    %% load files   
    [path1,file1,~] = fileparts(pathSimResults);
    addpath(path1);
    load(file1,'R');
        
    [path2,file2,~] = fileparts(pathExpData);
    addpath(path2);
    load(file2,'ExperimentalData');

    %% Helper names
    joints_ref = {'pelvis_tilt','pelvis_list','pelvis_rotation',...
        'hip_flexion','hip_adduction','hip_rotation',...
        'knee_angle','ankle_angle','subtalar_angle','tmt_angle','mtp_angle',...
        'lumbar_extension','lumbar_bending','lumbar_rotation',...
        'arm_flex','arm_add','arm_rot','elbow_flex'};
    joints_ref_r = {'pelvis_tilt','pelvis_list','pelvis_rotation',...
        'hip_flexion_r','hip_adduction_r','hip_rotation_r',...
        'knee_angle_r','ankle_angle_r','subtalar_angle_r','tmt_angle_r','mtp_angle_r',...
        'lumbar_extension','lumbar_bending','lumbar_rotation',...
        'arm_flex_r','arm_add_r','arm_rot_r','elbow_flex_r'};
    % since joint names are not consistent between measurement data and
    % simulation results :(
    joints_tit = {'Pelvis tilt','Pelvis list','Pelvis rotation','Pelvis tx',...
        'Pelvis ty','Pelvis tz','Hip flexion L','Hip adduction L',...
        'Hip rotation L','Hip flexion R','Hip adduction R','Hip rotation R',...
        'Knee L','Knee R','Ankle L','Ankle R',...
        'Subtalar L','Subtalar R','TMT L','TMT R','MTP L','MTP R',...
        'Lumbar extension','Lumbar bending','Lumbar rotation',...
        'Arm flexion L','Arm adduction L','Arm rotation L',...
        'Arm flexion R','Arm adduction R','Arm rotation R',...
        'Elbow flexion L','Elbow flexion R'};
    
        
         
       
    %% compare
    idx_Qs = [1,2,3,10,11,12,14,16,18,20,22,23,24,25,29,30,31,33];
       
    %ccc.joints =-ones(length(idx_Qs),1);
    ccc.kinematics.atShift0 = -ones(1,length(idx_Qs));      % -1 default, to see if field wasn't calculated
    ccc.kinematics.max = -ones(1,length(idx_Qs));           % -1 default, to see if field wasn't calculated
    ccc.kinematics.shiftAtMax = -ones(1,length(idx_Qs));    % -1 default, to see if field wasn't calculated

    ccc.kinetics.atShift0 = -ones(1,length(idx_Qs));      % -1 default, to see if field wasn't calculated
    ccc.kinetics.max = -ones(1,length(idx_Qs));           % -1 default, to see if field wasn't calculated
    ccc.kinetics.shiftAtMax = -ones(1,length(idx_Qs));    % -1 default, to see if field wasn't calculated

    subject = 'subject1';
    Qref = ExperimentalData.Q;
    IDref = ExperimentalData.Torques;
    
    for i = 1:length(idx_Qs)
        idx_jref = strcmp(Qref.(subject).Qs.colheaders,joints_ref{i});
        idx_sim = strcmp(R.colheaders.joints,joints_ref_r{i});
        ccc.joints{i} = joints_tit{idx_Qs(i)};
        if sum(idx_jref) == 1 && sum(idx_sim) == 1
            mean = Qref.(subject).Qs.mean(:,idx_jref);
            %std = Qref.(subject).Qs.std(:,idx_jref);
            sim = R.Qs(:,idx_sim);
            m = size(mean);
            s = size(sim);
            if(m == s)
            %rmse(i,:) = rsquare(mean,sim);
            %[c,lags] = xcorr(mean,sim,50);
            [c,lags] = xcorr(mean,sim,50,'coeff');     % normalized
            lags_i = find(lags == 0);
            ccc.kinematics.atShift0(1,i) = c(lags_i);   % ccc at shift = 0
            ccc.kinematics.max(1,i) = max(c);      % max ccc
            c_i = find(c == max(c));
            ccc.kinematics.shiftAtMax(1,i) = lags(c_i);   % shift of max ccc
            ccc.kinematics.shift(:,i) = lags';
            ccc.kinematics.ccc(:,i) = c';
            end

            mean = IDref.(subject).mean(:,idx_jref);
            %std = IDref.(subject).std(:,idx_jref);
            sim = R.Tid(:,idx_sim);
            m = size(mean);
            s = size(sim);
            if(m == s)
            %rmse(i,:) = rsquare(mean,sim);
            %[c,lags] = xcorr(mean,sim,50);
            [c,lags] = xcorr(mean,sim,50,'coeff');     % normalized
            lags_i = find(lags == 0);
            ccc.kinetics.atShift0(1,i) = c(lags_i);   % ccc at shift = 0
            ccc.kinetics.max(1,i) = max(c);      % max ccc
            c_i = find(c == max(c));
            ccc.kinetics.shiftAtMax(1,i) = lags(c_i);   % shift of max ccc
            ccc.kinetics.shift(:,i) = lags';
            ccc.kinetics.ccc(:,i) = c';
            end
        end
    end
    
  
    
    
    
end