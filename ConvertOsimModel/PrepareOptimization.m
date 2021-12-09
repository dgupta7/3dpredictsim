% This script extracts information from an OpenSim model, and stores it in
% a format that is accessible for formulating the simulation.
%
% Author: Lars D'Hondt
% Date: 7/Dec/2021
%

clear
clc

%% Inputs

% OpenSim file name
OsimFileName = 'subject1_mtj.osim';
% model selection options: Rajagopal, Gait92, Gait92_mtj
ModelName = 'Gait92_mtj';

% Boolean to select if we have to run the muscle analysis
Bool_RunMA = 1; 


%% Path information
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath(fullfile(pathRepo,'\VariousFunctions'));
% addpath(fullfile(pathRepo,'\CasADiFunctions'));
addpath(fullfile(pathRepo,'\Polynomials'));
% addpath(fullfile(pathRepo,'\MetabolicEnergy'));
% AddCasadiPaths();

% Modelpath
ModelPath = fullfile(pathRepo,'OpenSimModel/subject1',OsimFileName);

% Folder to save the polynomials
PolyFolder = OsimFileName(1:end-5);
pathPoly = fullfile(pathRepo,'Polynomials',PolyFolder);
if ~isfolder(pathPoly)
    mkdir(pathPoly);
end
pathMusc = fullfile(pathRepo,'MuscleModel',PolyFolder);
if ~isfolder(pathMusc)
    mkdir(pathMusc);
end

%% Fit polynomial functions to length and moment arms of muscles and ligaments
FitPolynomials(pathRepo,ModelName,ModelPath,PolyFolder,Bool_RunMA)
getPlantarFasciaGeometry(pathRepo,Modelpath,PolyFolder)

%% Get muscle parameters
% muscles are hard coded
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

MTparameters = getMTparameters(ModelPath,muscleNames);

save(fullfile(pathMusc,'MTparameters.mat'),'MTparameters');


















