%% hyperelastic plantar fascia model
% parameters obtained from: Kitaoka, HB; Luo, ZP; Growney, ES; Berglund, LJ;
% An, KN: Material properties of the plantar aponeurosis. Foot Ankle Int 15(10): 557 – 560, 1994.


lz = linspace(1,1.3,1000);
ls = 0.17;
l = lz*ls;


%Mooney-Rivlin Calculator
%Parameters are defined
c10=-222.1;
c01=290.97;
c20=-1.1257;
c11=4.7267;
c02=79.602;
%Paramters defined as variables, stretch defined as variable x
syms x %c10 c01 c20 c11 c02
% First and Second Invariants Defined
y=((x^2)+(2/x));
z=((2*x)+(x^-2));
%Strain Energy Defined
f1=((c10*(y-3))+(c01*(z-3))+(c11*(y-3)*(z-3))+((c20*((y-3)^2)))+((c02*((z-3)^2))));
%Derivative Taken
sigma=x*diff(f1);
%sigma@x= sigma(x) for any x value
f_sigma = matlabFunction(sigma);

sg = feval(f_sigma,lz);
sg_e = sg./lz;

c10 = -222.1;
c01 = 290.97;
c20 = -1.1257;
c11 = 4.7267;
c02 = 79.602;


for i=1:length(lz)
    I1 = lz(i)^2 + 2*lz(i)^(-1);
    I2 = lz(i)^(-2) + 2*lz(i);

    W(i) = c10*(I1-3) + c01*(I2-3) + c20*(I1-3)^2 + c11*(I1-3)*(I2-3) + c02*(I2-3)^2;
    
%     sig(i) = c10*(2*lz(i)-2/lz(i)^2) + c01*(2-2/lz(i)^3) + 4*c20*(lz(i)^3+1-2/lz(i)^3) +...
%         2*c11*(3*lz(i)^2-1+1/lz(i)-1/lz(i)^3-2/lz(i)^4) + 4*c02*(2*lz(i)-1/lz(i)^2-1/lz(i)^5);

%     sig1(i) = c10 + 2*c20*(I1-3) + c11*(I2-3);
%     sig2(i) = c01 + c11*(I1-3) + 2*c02*(I2-3);

    
end





dW = (W(2:end) - W(1:end-1))./(lz(2:end) - lz(1:end-1));
dW = [0,dW];
s = dW; % engineering stress
T = lz.*dW; % Cauchy stress


% F = s * 290;
F = T * 290;


dl = l-ls;
k = F./dl;

%%

figure
hold on
plot(lz,sg)
plot(lz,sg_e,':')
plot(lz,s,'--')
plot(lz,T,'-.')

%%
% figure
% plot(dl,F)
% xlabel('l/ls')
% ylabel('F')

% figure
% plot(dl,k)
% xlabel('l/ls')
% ylabel('k')

figure
plot(lz,F)
hold on
tau = 0.0125;
F_appr = 11e3*(1-exp((1-lz)*0.14/tau));
plot(lz,F_appr)
xlabel('l/ls')
ylabel('F')

% figure
% plot(lz,W)
% xlabel('l/ls')
% ylabel('W')

%%

a1 = -488737.9;
a2 = 2648898.5;
a3 = -5736967.6;
a4 = 6206986.7;
a5 = -3354935.1;
a6 = 724755.5;

sgm = a1*lz.^5 + a2*lz.^4 + a3*lz.^3 + a4*lz.^2 + a5*lz + a6;

Pr = 0.4; % Poisson ratio
A0 = 290; % initial section (mm^2)
A = A0*(1-Pr*(lz-1)).^2;

F = sgm.*A;

% figure
% plot(lz,F)
% hold on
% plot(lz,sgm*A0)
% % plot(lz,7e5*(l-ls))
% xlabel('l/ls')
% ylabel('F')

E = sgm./(lz-1);

% figure
% plot(lz(2:end),E(2:end))
% xlabel('l/ls')
% ylabel('E')


%% Plantar fascia stiffness model

% data from DOI: 10.3390/app11041517
% method from https://doi.org/10.1016/j.jbiomech.2017.10.037


Es = [11.82, 10.52, 10.70, 14.02, 13.85, 9.72, 8.78, 14.43, 8.42, 9.87, 18.7, 10.28, 6.75, 9.99, 13.17];
sigma = mean(Es)*0.08; % mean stress at 8% strain
e_0 = 0.055;
sigma_0 = 0.40;

E = (sigma)/(0.08-e_0);

% initial dimentions
ls = 0.17;
A0 = 290;

% parameter identification
mu = -ls*e_0;
F0 = sigma_0*A0;
k = E*A0/ls;
std = sqrt(2*pi)*F0/k;


x = linspace(-2,20,100)'*1e-3;

z = (x+mu)/(sqrt(2)*std);


d_erf = @(t) exp(-t.^2);
% d_erf_1 = exp(-z.^2);

for i=1:length(x)
   erf(i,1) = 2/sqrt(pi)* integral(d_erf,0,z(i));
   
   zi = linspace(0,z(i),200)';
   d_erf_i = exp(-zi.^2);
   erf(i,2) = 2/sqrt(pi)* trapz(d_erf_i)*z(i)/200;
end

% figure
% plot(z,erf)


F1 = k*std/sqrt(2*pi)*exp(-(x+mu).^2/(2*std^2));
F2 = k*(x+mu)/2 .*(erf+1);
F = F1+F2;

% figure
% plot(x,F1)
% hold on
% plot(x,F2)
% plot(x,F)

figure
plot((x/ls),F/A0)
xlim([0,0.1])
ylim([0,5])
