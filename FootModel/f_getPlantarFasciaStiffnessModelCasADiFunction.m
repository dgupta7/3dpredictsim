function f_PF_stiffness = f_getPlantarFasciaStiffnessModelCasADiFunction(modelType,varargin)
% This function returns a casADi function object that describes the tensile
% force of the plantar fascia in function of its length.
%
% modelType contains a string specifying the model to be used. Valid
% entries can be seen in section 3).
%
% It is possible to overwrite default parameter values by adding their
% name (string) and desired value (double) as additional inputs. Do note
% that most model types use different parameters.
%
% Author: Lars D'Hondt (March 2021)


% variable
import casadi.*
l = SX.sym('l');

%% 1) Default parameters
E = 350; % Young's modulus (N/mm^2)
% A0 = 290; % initial cross-section (mm^2)
A0 = 49.7; % initial cross-section (mm^2)
nu = 0.4; % Poisson ratio
ls = 0.135; % slack length (m)
k = E*A0/ls; % spring constant (N/m)
dl_0 = ls/100*1; % toe-in length

%% 2) Check for possible non-default parameter inputs
if length(varargin) >= 2 && mod(length(varargin),2) == 0
   for i=1:2:length(varargin)
       if strcmp(varargin{i},'E')
           E = varargin{i+1};
       elseif strcmp(varargin{i},'A0')
           A0 = varargin{i+1};
       elseif strcmp(varargin{i},'ls')
           ls = varargin{i+1};
       end
       % might need to recalculate those:
       k = E*A0/ls; % spring constant (N/m)
       dl_0 = ls/100*1; % toe-in length
       if strcmp(varargin{i},'k1')
           k = varargin{i+1};
       elseif strcmp(varargin{i},'dl_0')
           dl_0 = varargin{i+1};
       end
   end
end

%% 3) Calculate force 
% intermediate variables
cf = (tanh((l-1.0005*ls)/ls*1e4)+1)/2; 
dl_1 = l-ls; % elongation
dl = dl_1*cf; % positive elongation
lambda = l/ls; % stretch ratio
A = A0*lambda^(-nu*2); % actual cross-section

% stiffness models
if strcmp(modelType,'linear')
    F_PF = k*dl;
    
elseif strcmp(modelType,'tanh')
    F_PF = k*(dl - dl_0*tanh(dl/dl_0));
    
elseif strcmp(modelType,'Gefen2001')
    % 5th order polynomial approximation
    % https://doi.org/10.1016/S0021-9290(01)00242-1
    a1 = -488737.9;
    a2 = 2648898.5;
    a3 = -5736967.6;
    a4 = 6206986.7;
    a5 = -3354935.1;
    a6 = 724755.5;
    sigma = a1*lambda^5 + a2*lambda^4 + a3*lambda^3 + a4*lambda^2 + a5*lambda + a6 -0.100; % stress
    F_PF = sigma*A; % correction term, to make F=0 for l=ls
    
elseif strcmp(modelType,'Cheng2008')
    % Mooney-Rivlin model with 5 parameters (2nd order)
    % DOI: 10.3113/FAI.2008.0845
    % https://digitalcommons.unf.edu/etd/740
    c10 = -222.1;
    c01 = 290.97;
    c20 = -1.1257;
    c11 = 4.7267;
    c02 = 79.602;
    lambda = l/ls; % stretch ratio
    I1 = ((lambda^2)+(2/lambda)); % 1st strain invariant
    I2 = ((2*lambda)+(lambda^-2)); % 2nd strain invariant
    W = c10*(I1-3) + c01*(I2-3) + c11*(I1-3)*(I2-3) + c20*(I1-3)^2 + ...
        c02*(I2-3)^2; % strain energy density function
        % sigma = lambda * jac(W,lambda)
        % sigma_e = sigma/lambda, assuming uniaxial tension
        % combining both steps gives:
    sigma_e = ls*jacobian(W,l); % engineering stress
    F_PF = sigma_e*A0;
    
elseif strcmp(modelType,'Barrett2018') || strcmp(modelType,'toein_gaussian')
    % data from DOI 10.1007/s00276-011-0873-z
    e_0 = 0.03; % nominal strain at 0 stress when extrapolating linear region (-)
    sigma_0 = 0.5; % nominal strain at e_0 (MPa)
    E = 72;
    % method from https://doi.org/10.1016/j.jbiomech.2017.10.037
    mu = -ls*e_0;
    F0 = sigma_0*A0;
    k = E*A0/ls;
    std = sqrt(2*pi)*F0/k;
    x = l - ls;
    z = (x+mu)/(sqrt(2)*std);
    zi = linspace(0,z,500)';
    d_erf = exp(-zi.^2);
    erf = 2/sqrt(pi)* sum(d_erf(2:end)*z)/500;
    F_PF = k*std/sqrt(2*pi)*exp(-(x+mu).^2/(2*std^2)) + k*(x+mu)/2 .*(erf+1);
    
elseif strcmp(modelType,'Natali2010')
    % https://doi.org/10.3109/03008200903389127
%     mu = 14.449; % (MPa)
    mu = 0;
    k = 254.02; % (MPa)
    alpha = 10.397; % (-)
    sigma = mu*(lambda^2 - 1/lambda) + k/(2*alpha) *(exp(alpha*(lambda^2-1))-1)*lambda^2; % Cauchy stress
    F_PF = sigma*A;
    
elseif strcmp(modelType,'Ker1987')
    % see ligament_torques_Ker87.m
%     F_PF = -3578055333 + 21011527142*lambda^1 + -51373584416*lambda^2 + 66942934105*lambda^3 + ...
%         -49031631482*lambda^4 + 19139494724*lambda^5 + -3110684740.1*lambda^6; %sf=1, corrf

%     F_PF = -17217144166 + 98866378186.7*lambda^1 + -236477798399*lambda^2 + 301574344494*lambda^3 +...
%         -216263086165*lambda^4 + 82685305006.2*lambda^5 +
%         -13167998957.6*lambda^6; %sf=0.9

    F_PF = -15561008725 + 89385254015.1*lambda^1 + -213868002470*lambda^2 + 272825909452*lambda^3 +...
        -195707300417*lambda^4 + 74848663484.9*lambda^5 + -11923515339.9*lambda^6; %sf=0.9, F*sf^2

else
    % plantar fascia fully released
    F_PF = 0;
    warning('No valid function to describe the plantar fascia stiffness model, using 0 instead.')
end  

F = F_PF*cf; % F >= 0


%% 4) Build function
f_PF_stiffness = Function('f_PF_stiffness',{l},{F},{'l'},{'F'});

















end