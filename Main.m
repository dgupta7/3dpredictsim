%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script serves as the main file for the changes made in the branch 
% Thesis_Lars. It allows to specify the settings, solve, and post-process  
% a single gait simulation.
%
% Alternatively, the specified settings can be saved and added to a batch.
% To run the batch, run \BatchRunQueue.m
%
% To plot multiple results on the same figure, use \Plots\make_any_plot.m
%
% Run \ConvertOsimModel\PrepareOptimization.m 
%
% Author: Lars D'Hondt
% Date: December 2021
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
close all
clc

%% Paths
[pathRepo,~,~] = fileparts(mfilename('fullpath'));
addpath([pathRepo '/OCP']);
addpath([pathRepo '/VariousFunctions']);
addpath([pathRepo '/Plots']);
addpath([pathRepo '/CasADiFunctions']);
addpath([pathRepo '/Musclemodel']);
addpath([pathRepo '/Polynomials']);
addpath([pathRepo '/Debug']);
addpath([pathRepo '/FootModel']);
addpath([pathRepo '/RunSim']);
% AddCasadiPaths();

%% General settings
%-------------------------------------------------------------------------%
% Full body gait simulation
run_simulation = 1;         % run solver
post_process_results = 0;   % postproces
add_to_batch_queue = 0;     % save settings to run later

% settings for optimization
S.v_tgt     = 1.33;     % average speed
S.N         = 50;       % number of mesh intervals
S.NThreads  = 6;        % number of threads for parallel computing
S.max_iter  = 10;       % maximum number of iterations (comment -> 10000)
% S.tol_ipopt = 3;        % stopping criterion: < 10^(-...) 

% output folder
S.ResultsFolder = 'debug'; % subfolder of \Results where the result will be saved
% S.suffixCasName = 'v3';     % suffix for name of folder with casadifunctions
% S.suffixName = 'v4';        % suffix for name of file with results


%% Foot model
%-------------------------------------------------------------------------%
% General
S.Foot.Model = 'mtj';
   % 'mtp': foot with mtp joint
   % 'mtj': foot with mtp and midtarsal joint
S.Foot.Scaling = 'custom'; % default, custom, personalised

% Contact spheres
S.Foot.contactStiffnessFactor = 10;  % 1 or 10, 10: contact spheres are 10x stiffer
S.Foot.contactSphereOffsetY = 1;    % contact spheres are offset in y-direction to match static trial IK
S.Foot.contactSphereOffset1X = 0;   % heel contact sphere offset in x-direction

%% metatarsophalangeal (mtp) joint
S.Foot.mtp_actuator = 0;    % use an ideal torque actuator
S.Foot.mtp_muscles = 1;     % extrinsic toe flexors and extensors act on mtp joint
S.Foot.kMTP = 1;            % additional stiffness of the joint (Nm/rad)
% S.Foot.dMTP = 0.3;          % additional damping of the joint (Nms/rad)

%% midtarsal joint 
% (only used if Model = mtj)
S.Foot.mtj_muscles = 1;  % joint interacts with- extrinsic foot muscles
% lumped ligaments (long, short planter ligament, etc)
S.Foot.MT_li_nonl = 1;       % 1: nonlinear torque-angle characteristic
S.Foot.mtj_stiffness = 'Gefen2002';

S.Foot.kMT_li = 200;        % angular stiffness in case of linear
S.Foot.kMT_li2 = 10;        % angular stiffness in case of signed linear
% S.dMT = 1;                % (Nms/rad) damping

% plantar fascia
S.Foot.PF_stiffness = 'Gefen2002'; % 'none''linear''Gefen2002''Cheng2008''Natali2010''Song2011'
S.Foot.PF_sf = 1;   
S.Foot.PF_slack_length = 0.15; % (m) slack length

% Plantar Intrinsic Muscles represented by and ideal force actuator
S.Foot.PIM = 0;     % include PIM actuator
S.W.PIM = 10^3;     % weight on the excitations for cost function
S.W.P_PIM = 10^4;   % weight on the net Work for cost function


%% Initial guess
%-------------------------------------------------------------------------%

% initial guess based on simulations without exoskeletons
    % initial guess identifier (1: quasi random, 2: data-based)
S.IGsel         = 2;
    % initial guess mode identifier (1 walk, 2 run, 3 prev.solution, 4 solution from /IG/Data folder)
S.IGmodeID      = 1;

if S.IGmodeID == 4
    S.savename_ig   = 'NoExo';
elseif S.IGmodeID == 3
    S.ResultsF_ig   = '';
    S.savename_ig   = 'Fal_s1_bCst_ig21';
end


%% run simulation
PredSim(S,run_simulation,post_process_results,add_to_batch_queue);














