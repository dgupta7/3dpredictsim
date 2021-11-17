clear
clc

% soleus_r
lMT_mtp = 0.281;
lMT_mtj = 0.294;

lM0_mtp = 0.048170551124502803;
lTs_mtp = 0.24085275562251399;

sf = lMT_mtj/lMT_mtp;
% sf = sf*0.99;

lM0_mtj = lM0_mtp*sf;
lTs_mtj = lTs_mtp*sf;

disp(['lM0: ' num2str(lM0_mtj,5)])
disp(['lTs: ' num2str(lTs_mtj,5)])

% lat_gas_r

% med_gas_r

% tib_post_r


% tib_ant_r

% per_brev_r

% per_long_r

% per_tert_r

% flex_dig_r

% flex_hal_r

% ext_dig_r

% ext_hal_r


