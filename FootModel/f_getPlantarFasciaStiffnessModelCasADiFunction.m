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
A0 = 290; % initial cross-section (mm^2)
ls = 0.17; % slack length (m)
k1 = E*A0/ls; % spring constant (N/m)
k2 = 5e7; % spring constant (N/m^2)
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
       k1 = E*A0/ls; % spring constant (N/m)
       dl_0 = ls/100*1; % toe-in length
       if strcmp(varargin{i},'k1')
           k1 = varargin{i+1};
       elseif strcmp(varargin{i},'k2')
           k2 = varargin{i+1};
       elseif strcmp(varargin{i},'dl_0')
           dl_0 = varargin{i+1};
       end
   end
end

%% 3) Calculate force 
if strcmp(modelType,'linear')
    F_PF = k1*(l-ls);
    
elseif strcmp(modelType,'hypoelastic_tanh')
    F_PF = k1*((l-ls) - dl_0*tanh((l-ls)/dl_0));
    
elseif strcmp(modelType,'hypoelastic_sqr')
    F_PF = k2*(l-ls)^2;
    
elseif strcmp(modelType,'hypoelastic_poly5')
    % 5th order polynomial approximation
    % https://doi.org/10.1016/S0021-9290(01)00242-1
    a1 = -488737.9;
    a2 = 2648898.5;
    a3 = -5736967.6;
    a4 = 6206986.7;
    a5 = -3354935.1;
    a6 = 724755.5;
    Pr = 0.4; % Poisson ratio
    lambda = l/ls; % stretch ratio
    sigma = a1*lambda^5 + a2*lambda^4 + a3*lambda^3 + a4*lambda^2 + a5*lambda + a6; % stress
    A = A0*(1-Pr*(lambda-1))^2; % actual cross-section
    F_PF = sigma*A;
    
elseif strcmp(modelType,'hyperelastic_MR5')
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
    
end

F = F_PF.*( tanh(F_PF)+1 )/2; % F >= 0


%% 4) Build function
f_PF_stiffness = Function('f_PF_stiffness',{l},{F},{'l'},{'F'});

















end