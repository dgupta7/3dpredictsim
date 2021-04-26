

ExtPoly = '_mtp';
subject = R.S.subject;

pathmusclemodel = fullfile(pathRepo,'MuscleModel',subject);
addpath(genpath(pathmusclemodel));

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

musi = 1:length(muscleNames(1:end-3));

load([pathmusclemodel,'/MTparameters_',subject, ExtPoly, '.mat']);
MTparameters_m = [MTparameters(:,musi),MTparameters(:,musi)];

NMuscle = 92;
aTendon = 35*ones(NMuscle,1);
IndexCalf = [32 33 34 78 79 80];    % adjust stiffness of the calf muscles
aTendon(IndexCalf) = 35;
shift = getShift(aTendon);

iSol = find(strcmp(R.colheaders.muscles,'soleus_r'));
m = iSol;

for i=1:100
    FTtilde = R.FTtilde(i,:);
    dFTtilde = R.dFTtilde(i,:);
    
    lMT = R.lMT(i,:);
    vMT = R.vMT(i,:);
    
    [vM0,vMtilde0,vT0] = FiberVelocity_TendonForce_tendon(FTtilde(m),...
            dFTtilde(m),MTparameters_m(:,m),lMT(m),vMT(m),aTendon(m),shift(m));

        vM(i) = vM0;
        vMtilde(i) = vMtilde0;
        vT(i) = vT0;
        P_ach(i) = -vT(i)*R.FT(i,iSol);
end



%%
figure
plot(P_ach)






