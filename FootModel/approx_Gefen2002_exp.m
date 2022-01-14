

lambda1 = linspace(1,1.15,500)';
lambda2 = linspace(1,1.2,500)';
lambda3 = linspace(0.95,1.27,500)';

figure
p1=plot(lambda3,getLigamentGefen2002(lambda3),'DisplayName','5e order polynomial');
hold on
xlabel('\lambda (-)','Interpreter','tex')
ylabel('\sigma (N/mm^2)','Interpreter','tex')
% xlim([0.99,lambda3(end)])
ylim([-10,70])

%%
modelfun_G_exp = @(coeff,x) coeff(1)*(exp( coeff(2)*(x-coeff(3)) ) - 1);
coeff_0 = [6.7,9,1];
sigma_fit = modelfun_G_exp(coeff_0,lambda2);
% plot(lambda2,sigma_fit)

%%



x_fit = lambda1;
y_fit = getLigamentGefen2002(lambda1);

mdl = fitnlm(x_fit,y_fit,modelfun_G_exp,coeff_0);
coeff_sol = table2array(mdl.Coefficients(:,1));
f_getMtjLigamentMoment = @(lm) modelfun_G_exp(coeff_sol,lm);

p2=plot(lambda3,f_getMtjLigamentMoment(lambda3),'DisplayName','c_1(exp(c_2[\lambda - c_3]) - 1)');
xline(1,'-k')
xline(1.15,'-k','fitting range','LabelHorizontalAlignment','left','LabelOrientation','horizontal')
legend([p1,p2],'Location','southeast','Interpreter','tex')
title('Ligament stiffness (Gefen, 2002)')

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
%     sigma_f(lambda_f<1) = 0;
end