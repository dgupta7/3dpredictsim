
clear
close all
clc

pathmain = pwd;
[pathRepo,~,~]  = fileparts(pathmain);


%%
Modelpath = [pathRepo '\OpenSimModel\subject1\DetailedFootModel_pinMTJ_sd.osim'];
PolyFolder = [pathRepo '\Polynomials\Fal_s1_mtj_sd'];

ligament_names = {'LongPlantar1_r','LongPlantar2_r','LongPlantar3_r','LongPlantar4_r',... % long plantar ligament
    'CalcaneoCuboidPlantar1Mus_r','CalcaneoCuboidPlantar2Mus_r',... % short plantar ligament
    'CalcaneoNavicularPlantar1Mus_r','CalcaneoNavicularPlantar2Mus_r','CalcaneoNavicularPlantar3Mus_r',... % spring ligament
    'CalcaneoCuboidDorsalMus_r','CalcaneoNavicularBifurcateMus_r','CalcaneoCuboidBifurcateMus_r'}; % other


p = length(ligament_names);
n = 500;

q_mtj = linspace(-25,15,n)*pi/180;


%% Initialise model
import org.opensim.modeling.*;
model = Model(Modelpath);
s = model.initSystem;


for i=0:model.getCoordinateSet().getSize()-1
    if strcmp(char(model.getCoordinateSet().get(i)),'mtj_angle_r')
        idx_mtj = i;
        break
    end
end


%% Evaluate ligament moment for dummy motion
% Set state vector to 0
state_vars = model.getStateVariableValues(s);
state_vars.setToZero();
model.setStateVariableValues(s,state_vars);
model.realizePosition(s);
% Initialise matrices for results
l_0 = zeros(1,p);
pcsa_0 = zeros(1,p);
l_lig = zeros(n,p);
d_lig = zeros(n,p);
F_lig = zeros(n,p);

for j=1:p
    lig = Ligament.safeDownCast(model.getForceSet().get(ligament_names{j}));
    l_0(1,j) = lig.getRestingLength();
    pcsa_0(1,j) = lig.get_pcsa_force();
    Fl_c = lig.get_force_length_curve();

    for i=1:n
        % Set each coordinate value
        state_vars.set(idx_mtj*2,q_mtj(i));
        model.setStateVariableValues(s,state_vars);
        model.realizePosition(s);
        % Get length
        l_lig(i,j) = lig.getLength(s);
        % Get moment arm
        d_lig(i,j) = lig.computeMomentArm(s,model.getCoordinateSet().get(idx_mtj));
        
        % Get force
        l_ij = Vector(1,l_lig(i,j)/l_0(1,j));
        F_ij = Fl_c.calcValue(l_ij);
        F_lig(i,j) = F_ij*pcsa_0(1,j);

    end

end

M_lig = d_lig.*F_lig;
M_mtj = sum(M_lig,2);


%%
CsV = hsv(p);
mrk = {'-','--',':','-.'};

fg1=figure;
subplot(3,1,1)
plot(q_mtj*180/pi,M_mtj)
xlabel('Mtj angle (°)')
ylabel('Mtj logament moment (Nm)')

for i=1:p
    subplot(3,1,2)
    plot(q_mtj*180/pi,M_lig(:,i),'Color',CsV(i,:),'LineStyle',mrk{rem(i,4)+1})
    hold on
    xlabel('Mtj angle (°)')
    ylabel('Mtj logament moment (Nm)')
    
    subplot(3,1,3)
    plot(l_lig(:,i)/l_0(i),F_lig(:,i),'Color',CsV(i,:),'LineStyle',mrk{rem(i,4)+1})
    hold on
    xlim([0.99,1.1])
    xlabel('\lambda (-)')
    ylabel('Tension (N)')
end
legend(ligament_names,'Interpreter','none','Location','northeastoutside')

%%
fg2=figure;
for i=1:p
    e_lig = l_lig(:,i)/l_0(i)-1;
    e_lig_all(:,i) = e_lig;
    e_idx = find(e_lig>=0 & e_lig<=0.1);
    e_lig = e_lig(e_idx);
%     disp(length(e_lig))
%     p1 = polyfit(e_lig,F_lig(e_idx,i),1);
%     EA(i) = p1(1);

    plot(e_lig+1,F_lig(e_idx,i)/pcsa_0(i))
    hold on
    tmp = F_lig(e_idx,i);
    if isempty(tmp)
        tmp(1) = 0;
    end
    F_rel_max(i) = tmp(end)/pcsa_0(i);

    idx_tmp = find(e_lig_all(:,i)<0.2,1,"last");
    if isempty(idx_tmp)
        idx_tmp(1) = n;
    end
    idx_max(i) = idx_tmp;
end
q_max = q_mtj(min(idx_max))*180/pi;
E_mean_M = max(F_rel_max)/0.1;


lambda = linspace(1,1.1,500)';
sigma = getLigamentGefen2002(lambda);

nu = 0.4;

sigma_eng = sigma.*lambda.^(-nu*2);

fg3=figure;
plot(lambda,sigma)
hold on
plot(lambda,sigma_eng)

pol =  polyfit(lambda-1,sigma_eng,1);
E_mean_G = pol(1);

sigma_pol = polyval(pol,lambda-1);
plot(lambda,sigma_pol,'--')

E_mean_G1 = (lambda-1)\sigma_eng;
plot(lambda,E_mean_G1*(lambda-1),'--')
plot([1,1.1],[0,0.1*E_mean_G1],'--')


sf = E_mean_M/E_mean_G1;
figure(fg2)
hold on
plot([1,1.1],[0,0.1*E_mean_M],'--')
plot([1,1.1],[0,0.1*E_mean_G1*sf],'--')
plot(lambda,sigma_eng*sf,'--k')

%%
F_lig_G = zeros(size(F_lig));

for j=1:p
    lambda_j = l_lig(:,j)/l_0(1,j);
    sigma_j = getLigamentGefen2002(lambda_j);
    F_lig_G(:,j) = sigma_j*sf*pcsa_0(j);
end

M_lig_G = d_lig.*F_lig_G;
M_mtj_G = sum(M_lig_G,2);


fg4=figure;
subplot(3,1,1)
plot(q_mtj*180/pi,M_mtj_G)
hold on
plot(q_mtj*180/pi,M_mtj,'--k')
xlabel('Mtj angle (°)')
ylabel('Mtj logament moment (Nm)')

for i=1:p
    subplot(3,1,2)
    plot(q_mtj*180/pi,M_lig_G(:,i),'Color',CsV(i,:),'LineStyle',mrk{rem(i,4)+1})
    hold on
    xlabel('Mtj angle (°)')
    ylabel('Mtj logament moment (Nm)')
    
    subplot(3,1,3)
    plot(l_lig(:,i)/l_0(i),F_lig_G(:,i),'Color',CsV(i,:),'LineStyle',mrk{rem(i,4)+1})
    hold on
    xlim([0.99,1.1])
    xlabel('\lambda (-)')
    ylabel('Tension (N)')
end
legend(ligament_names,'Interpreter','none','Location','northeastoutside')

%%

fg5=figure;
plot(q_mtj*180/pi,M_mtj,'DisplayName','Malaquias et al. (2017)')
hold on
plot(q_mtj*180/pi,M_mtj_G,'DisplayName','Gefen (2002)')
legend('Location','best','Interpreter','none')


%%
% modelfun = @(coeff,x) coeff(1) + coeff(2)*exp((coeff(3)+x).*coeff(4)) + coeff(5).*x;
% 
% coeff_0 = [0,-0.3,0.2,10,-10];
% M_tmp = modelfun(coeff_0,q_mtj);
% 
% figure(fg5)
% plot(q_mtj*180/pi,M_tmp,'--','DisplayName','IG a + b*exp{(c+q)*d} + e*q')

%%
% idx_fit = find(q_mtj*180/pi<=q_max);
% x_fit = q_mtj(idx_fit);
% y_fit = M_mtj_G(idx_fit);
% 
% 
% mdl = fitnlm(x_fit,y_fit,modelfun,coeff_0);
% coeff_sol = table2array(mdl.Coefficients(:,1));
% f_getMtjLigamentMoment = @(q) modelfun(coeff_sol,q);


%%

% M_mtj_func = f_getMtjLigamentMoment(q_mtj);
% 
% figure(fg5)
% plot(q_mtj*180/pi,M_mtj_func,'--','DisplayName','a + b*exp{(c+q)*d} + e*q')

%%

% poly_GM = polyfit(q_mtj,M_mtj_G,5);
% figure(fg5)
% plot(q_mtj*180/pi,polyval(poly_GM,q_mtj),'--','DisplayName','5th order polynomial')

%%
% coeff = poly_GM(end:-1:1);
% codeString = ['M_li = ' num2str(coeff(1),5)];
% for ii=1:length(coeff)-1
%     codeString = [codeString ' + ' num2str(coeff(ii+1),5) '*q_mtj.^' num2str(ii)];
% end
% codeString = [codeString ';'];
% disp(codeString)

%%
% M_li = 0.081222 + -10.873*q_mtj.^1 + -307.95*q_mtj.^2 + -3487.8*q_mtj.^3 + -7202.7*q_mtj.^4 + 76515*q_mtj.^5;
% figure(fg5)
% plot(q_mtj*180/pi,M_li,'-.','DisplayName','5e-th order polynomial')

%%
import casadi.*
f_getMtjLigamentMoment = interpolant('f_getMtjLigamentMoment','bspline',{q_mtj},M_mtj_G);

figure(fg5)
M_mtj_cas = f_getMtjLigamentMoment(q_mtj);
plot(q_mtj*180/pi,full(M_mtj_cas),'--','DisplayName','casadi function')

f_getMtjLigamentMoment.save((fullfile(PolyFolder,'f_getMtjLigamentMoment')));

%%
import casadi.*
f_test = Function.load((fullfile(PolyFolder,'f_getMtjLigamentMoment')));

figure(fg5)
M_mtj_cas_t = f_test(q_mtj);
plot(q_mtj*180/pi,full(M_mtj_cas_t),'--','DisplayName','casadi function 1')

%%

function sigma_f = getLigamentGefen2002(lambda_f)
    % Gefen2002
    a1 = -412640.5;
    a2 = 2235967.7;
    a3 = -4841544.8;
    a4 = 5236972.7;
    a5 = -2829945.7;
    a6 = 611190.6;
    sigma_f = a1*lambda_f.^5 + a2*lambda_f.^4 + a3*lambda_f.^3 + a4*lambda_f.^2 + a5*lambda_f + a6;
    sigma_f(lambda_f<1) = 0;
end
