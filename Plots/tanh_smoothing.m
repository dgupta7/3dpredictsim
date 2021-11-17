clear
close all
clc
[pathHere,~,~] = fileparts(mfilename('fullpath'));
[pathRepo,~,~] = fileparts(pathHere);
addpath([pathRepo '/VariousFunctions']);
AddCasadiPaths();


%%
vM = linspace(-0.25,0.27,1000);

x = linspace(-5,5,1000);

figure
subplot(2,1,1)
plot(x,tanh(x))
xlabel('x')
ylabel('y')
title('y = tanh(x)')
grid on
ylim([-1.2,1.2])

subplot(2,1,2)
plot(x,0.5-0.5*tanh(x))
xlabel('x')
ylabel('y')
title('y = 0.5 - 0.5 tanh(x)')
grid on
ylim([-0.2,1.2])

%%

import casadi.*
x_SX = SX.sym('x');
b_SX = SX.sym('b');
xn_SX = (0.5 - 0.5*tanh(b_SX*x_SX));
y_SX = xn_SX*x_SX;
g_SX = gradient(y_SX,x_SX);

g1 = Function('g1',{x_SX,b_SX},{y_SX,g_SX});


%%
x = linspace(-pi,pi,1000);

x_neg = x;
x_neg(x>0) = 0;
g_neg = x_neg./x;

for i=1:length(x)
    [y_tmp,g_tmp] = g1(x(i),1);
    y_b1(i) = full(y_tmp);
    g_b1(i) = full(g_tmp);

    [y_tmp,g_tmp] = g1(x(i),2);
    y_b2(i) = full(y_tmp);
    g_b2(i) = full(g_tmp);

    [y_tmp,g_tmp] = g1(x(i),10);
    y_b10(i) = full(y_tmp);
    g_b10(i) = full(g_tmp);

end

figure
subplot(1,2,1)
plot(x,x_neg,'DisplayName','if(x < 0) y = x; else y = 0')
hold on
grid on
plot(x,y_b1,'DisplayName','y = [0.5 - 0.5 tanh(x)] x')
plot(x,y_b2,'DisplayName','y = [0.5 - 0.5 tanh(2x)] x')
plot(x,y_b10,'DisplayName','y = [0.5 - 0.5 tanh(10x)] x')
legend('Location','best')
ylabel('y')
xlabel('x')

subplot(1,2,2)
plot(x,g_neg,'DisplayName','if(x < 0) y = x; else y = 0')
hold on
grid on
plot(x,g_b1)
plot(x,g_b2)
plot(x,g_b10)
ylabel('\nabla_xy')
xlabel('x')

%%

load(fullfile([pathRepo '\Results\Final\Fal_s1_bCst_tanh10_ig21_pp.mat']),'R');
iM = [47:92];

vMs = R.Muscle.vM(:,iM(1));
vMTs = R.vMtilde(:,iM(1));
for ii=2:length(iM)
    vMs = [vMs; R.Muscle.vM(:,iM(ii))];
    vMTs = [vMTs; R.vMtilde(:,iM(ii))];
end

x = vM;

x_neg = x;
x_neg(x>0) = 0;
g_neg = x_neg./x;

for i=1:length(x)
    [y_tmp,g_tmp] = g1(x(i),10);
    y_b10(i) = full(y_tmp);
    g_b10(i) = full(g_tmp);

    [y_tmp,g_tmp] = g1(x(i),50);
    y_b50(i) = full(y_tmp);
    g_b50(i) = full(g_tmp);

    [y_tmp,g_tmp] = g1(x(i),100);
    y_b100(i) = full(y_tmp);
    g_b100(i) = full(g_tmp);

end

figure
subplot(2,1,1)
histogram(vMs,100)
xlabel('vM')
title('Fibre velocity of all muscles')
x_l = get(gca,'XLim');
y_l = get(gca,'YLim');
ylim([y_l(1),1.2*y_l(2)])

subplot(2,1,2)
plot(x,x_neg,'DisplayName','exact')
hold on
grid on
plot(x,y_b10,'DisplayName','b = 10')
plot(x,y_b50,'DisplayName','b = 50')
plot(x,y_b100,'DisplayName','b = 100')
legend('Location','best')
ylabel('vM negative')
xlabel('vM')
xlim(x_l)


figure
subplot(2,1,1)
histogram(vMTs,100)
xlabel('vMtilde')
title('Relative fibre velocity of all muscles')
x_l = get(gca,'XLim');
y_l = get(gca,'YLim');
ylim([y_l(1),1.2*y_l(2)])

subplot(2,1,2)
plot(x,x_neg,'DisplayName','exact')
hold on
grid on
plot(x,y_b10,'DisplayName','b = 10')
plot(x,y_b50,'DisplayName','b = 50')
plot(x,y_b100,'DisplayName','b = 100')
legend('Location','best')
ylabel('vM negative')
xlabel('vMtilde')
xlim(x_l)


%%

figure
plot(x,x_neg,'DisplayName','exact')
hold on
grid on
plot(x,y_b10,'DisplayName','b = 10')
plot(x,y_b50,'DisplayName','b = 50')
plot(x,y_b100,'DisplayName','b = 100')
legend('Location','best')
ylabel('vM negative')
xlabel('vM')
xlim([-0.035,0.054])
title('Soleus')


%%

N = hist(vMs,x(1:10:end));
figure
plot(x(1:10:end),(y_b10(1:10:end)-x_neg(1:10:end)).*N,'DisplayName','b = 10')
hold on
plot(x(1:10:end),(y_b50(1:10:end)-x_neg(1:10:end)).*N,'DisplayName','b = 50')
plot(x(1:10:end),(y_b100(1:10:end)-x_neg(1:10:end)).*N,'DisplayName','b = 100')
legend('Location','best')
xlabel('vM')
title('Weighted error')















