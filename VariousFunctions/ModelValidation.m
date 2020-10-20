function [rmse,ccc] = ModelValidation(pathSimResults, pathExpData)
% 


    %% load files   
    [path1,file1,~] = fileparts(pathSimResults);
    addpath(path1);
    load(file1);
        
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
    
        
        LegName = file1;
       
    %% compare kinematics
       
    idx_Qs = [1,2,3,10,11,12,14,16,18,20,21,22,23,27,28,29,31];
    rmse = -ones(length(idx_Qs),2);     % -1 default, to see if field wasn't calculated
    ccc = -ones(length(idx_Qs),3);     % -1 default, to see if field wasn't calculated
    
    
    subject = 'subject1';
    Qref = ExperimentalData.Q;
    figure()
    for i = 1:length(idx_Qs)
        idx_jref = strcmp(Qref.(subject).Qs.colheaders,joints_ref{i});
        mean = Qref.(subject).Qs.mean(:,idx_jref);
%         std = Qref.(subject).Qs.std(:,idx_jref);
        sim = R.Qs(:,idx_Qs(i));
        m = size(mean);
        s = size(sim);
        if(m == s)
        rmse(i,:) = rsquare(mean,sim);
        
        [c,lags] = xcorr(mean,sim,50,'coeff');
        lags_i = find(lags == 0);
        ccc(i,1) = c(lags_i);   % ccc at shift = 0
        ccc(i,2) = max(c);      % max ccc
        c_i = find(c == max(c));
        ccc(i,3) = lags(c_i);   % shift of max ccc
        
        subplot(3,6,i)
        stem(lags,c,'DisplayName',LegName)
        hold on
        grid on
        title(joints_tit{idx_Qs(i)})
        else
            disp(joints_ref{i});
        end

    end
    
    
    
    
end