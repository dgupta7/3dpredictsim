function [S] = GetDefaultSettings(S)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
if ~isfield(S,'savename_ig')
    S.savename_ig= [];
end

if ~isfield(S,'ResultsF_ig')
    S.ResultsF_ig = [];
end

if  ~isfield(S,'NThreads') || isempty(S.NThreads)
    S.NThreads = 4;
end

if ~isfield(S,'mass') || isempty(S.mass)
    S.mass = 62;
end

if ~isfield(S,'subject') || isempty(S.subject)
    S.subject = 'subject1';
end

% quasi random initial guess
if ~isfield(S,'IG_PelvisY') || isempty(S.IG_PelvisY)
    S.IG_PelvisY = 0.9385;
end

% Kinematics Constraints - Default Settings
if isfield(S,'Constr')
    if ~isfield(S.Constr,'calcn')
        S.Constr.calcn = 0.09;  % by default at least 9cm distance between calcn
    end
    if ~isfield(S.Constr,'toes')
        S.Constr.toes = 0.1; % by default at least 10cm distance between toes
    end
    if ~isfield(S.Constr,'tibia')
        S.Constr.tibia = 0.11; % by default at least 11cm distance between toes
    end
else
    S.Constr.calcn = 0.09;  % by default at least 9cm distance between calcn
    S.Constr.toes = 0.1; % by default at least 10cm distance between toes
    S.Constr.tibia = 0.11; % by default at least 11cm distance between toes
end
end

