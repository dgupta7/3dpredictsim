% function res = f_WindLassSlackParameters(vars)
vars = [0.14,0.2,1/50];

k_PF = 7.1994e+05;

l_s = vars(1);
q_s = vars(2);

dl_0 = 1e-3;
k_li = vars(3)*k_PF;

Qs_mtp = [-30:15:30]*pi/180;

cWL = 0.03;

a = 0.08207;
b = 0.089638;
phi0 = 2.493499;
H0 = 0.027280;

L0 = sqrt(a^2 + b^2 - 2*a*b*cos(phi0));



for i=1:length(Qs_mtp)
    q_mtp = Qs_mtp(i);
    l_0(i) = (1-cWL*(q_mtp*180/pi)/20)*L0; % foot arch length
    h_0(i) = a*b/l_0(i)*sin(phi0);
    q_mt_0(i) = acos( (a^2 + b^2 - l_0(i)^2)/(2*a*b) ) - phi0;

    dl_PF(i) = l_0(i) - l_s;
    F_PFi = PF_stiffness(k_PF,dl_PF(i),dl_0);
    F_PF(i) = F_PFi*( tanh(F_PFi)+1 )/2;
    
    F_li(i) = k_li*(q_mt_0(i)-q_s)/h_0(i);
    

    res(i) = F_PF(i) + F_li(i) - 5;
end
res

function force = PF_stiffness(k,dl,dl_stiff)
   force = k*(dl - dl_stiff*tanh(dl/dl_stiff));
end
% end