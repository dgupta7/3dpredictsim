
L = linspace(0.17,0.19,50)';


A0 = 290; % initial cross-section (mm^2)
ls = 0.17; % slack length (m)
a1 = -488737.9;
a2 = 2648898.5;
a3 = -5736967.6;
a4 = 6206986.7;
a5 = -3354935.1;
a6 = 724755.5;
nu = 0.4; % Poisson ratio

lambda = L/ls; % stretch ratio
sigma = a1*lambda.^5 + a2*lambda.^4 + a3*lambda.^3 + a4*lambda.^2 + a5*lambda + a6; % stress
e = lambda-1;

A = A0*(lambda.^(-nu)).^2;

F_PF = sigma.*A -29;

figure
plot(e,sigma)
hold on
plot(e,e*200)



a1+a2+a3+a4+a5+a6
