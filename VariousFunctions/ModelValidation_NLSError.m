function [NLSE] = ModelValidation_NLSError(pathSimResults, pathExpData)
% returns a struct with information about the nonlinear least squares error
% of the simulation results and the mean measurement data, weighted with
% the inverse of the variance


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
        'knee_angle','ankle_angle','subtalar_angle','mtp_angle',...
        'lumbar_extension','lumbar_bending','lumbar_rotation',...
        'arm_flex','arm_add','arm_rot','elbow_flex'};
    joints_tit = {'Pelvis tilt','Pelvis list','Pelvis rotation','Pelvis tx',...
        'Pelvis ty','Pelvis tz','Hip flexion L','Hip adduction L',...
        'Hip rotation L','Hip flexion R','Hip adduction R','Hip rotation R',...
        'Knee L','Knee R','Ankle L','Ankle R',...
        'Subtalar L','Subtalar R','MTP L','MTP R',...
        'Lumbar extension','Lumbar bending','Lumbar rotation',...
        'Arm flexion L','Arm adduction L','Arm rotation L',...
        'Arm flexion R','Arm adduction R','Arm rotation R',...
        'Elbow flexion L','Elbow flexion R'};
    
        
         
       
    %% compare
    idx_Qs = [1,2,3,10,11,12,14,16,18,20,21,22,23,27,28,29,31];
       

    NLSE.kinematics = -ones(1,length(idx_Qs));      % -1 default, to see if field wasn't calculated
    
    NLSE.kinetics = -ones(1,length(idx_Qs));        % -1 default, to see if field wasn't calculated
    
    subject = 'subject1';
    Qref = ExperimentalData.Q;
    IDref = ExperimentalData.Torques;
    
    for i = 1:length(idx_Qs)
        idx_jref = strcmp(Qref.(subject).Qs.colheaders,joints_ref{i});
        NLSE.joints{i} = joints_tit{idx_Qs(i)};
        mean = Qref.(subject).Qs.mean(:,idx_jref);
        std = Qref.(subject).Qs.std(:,idx_jref);
        sim = R.Qs(:,idx_Qs(i));
        m = size(mean);
        s = size(sim);
        if(m == s)
        NLSE.kinematics(1,i) = sum(((sim-mean)./std).^2);
        end
        
        mean = IDref.(subject).mean(:,idx_jref);
        std = IDref.(subject).std(:,idx_jref);
        sim = R.Tid(:,idx_Qs(i));
        m = size(mean);
        s = size(sim);
        if(m == s)
        NLSE.kinetics(1,i) = sum(((sim-mean)./std).^2);
        end
        

    end
    
  
    
    
    
end