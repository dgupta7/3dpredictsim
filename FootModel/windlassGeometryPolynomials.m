% approximating windlass geometry with polynomials
clear
close all
clc

q_mt = linspace(-30,30,100)'*pi/180;

%% Geometry
% run \FootModel\solveFootmodelParameters.m to get these values
% foot arch
calcn2mtj = 0.08207;
mtj2mtpj = 0.089638;
beta0 = 2.4935;
% windlass mechanism
calcnPF2mtj = 0.06695;
mtj2mttPF = 0.091714;
phi0 = 2.1274;

%% Geometry Windlass mechanism
% based on midtarsal joint angle
phi = phi0 + q_mt; % top angle of WL triangle
l_PF_fa = sqrt(calcnPF2mtj^2 + mtj2mttPF^2 - 2*calcnPF2mtj*mtj2mttPF.*cos(phi)); % length of PF spanning arch
MA_PF = calcnPF2mtj*mtj2mttPF./l_PF_fa.*sin(phi); % moment arm of PF to mtj



%% approximating length of PF spanning foot arch

scs = get(0,'ScreenSize');
figure('Position',[1+scs(3)/2,scs(4)/2+20,scs(3)/2, scs(4)/2-100]);
subplot(221)
plot(q_mt*180/pi,l_PF_fa*1e3,'DisplayName','exact')
hold on
grid on
legend('Location','southeast');
xlabel('Midtarsal angle (°)')
ylabel('l_{PF,fa} (mm)')
title('PF length spanning foot arch')
axis tight

%%
order = 2;
stop = 0;

while ~stop
    [l_fa,~,coeff1] = f_fitPolynomialFunction(q_mt,l_PF_fa,order);
    l = l_fa(q_mt);
    subplot(221)
    p1=plot(q_mt*180/pi,l*1e3,'--','DisplayName',[num2str(order) '^e order polynomial']);

    rel_err1 = abs((l-l_PF_fa)./l_PF_fa);
    err1 = max(rel_err1);
    subplot(223)
    semilogy(q_mt*180/pi,rel_err1,'Color',p1.Color)
    hold on

%     disp(err1)

    if err1 <= 1e-5 || order >= 5
        stop = 1;
    else
        order = order+1;
    end
end

subplot(223)
hold on
grid on
xlabel('Midtarsal angle (°)')
ylabel('log_{10}\deltal_{PF,fa} (-)')
title('Relative error of the approximation')
axis tight

codeString1 = ['l_PF_fa = ' num2str(coeff1(1),7)];
for ii=1:length(coeff1)-1
    codeString1 = [codeString1 ' + ' num2str(coeff1(ii+1),7) '.*q_mt.^' num2str(ii)];
end
codeString1 = [codeString1 ';'];
disp(codeString1)

%% approximating moment arm of PF around midtarsal joint

subplot(222)
plot(q_mt*180/pi,MA_PF*1e3,'DisplayName','exact')
hold on
grid on
legend('Location','northeast');
xlabel('Midtarsal angle (°)')
ylabel('MA_{PF} (mm)')
title('Moment arm PF')
axis tight

order = 2;
stop = 0;

while ~stop
    [MA,~,coeff2] = f_fitPolynomialFunction(q_mt,MA_PF,order);
    h = MA(q_mt);
    subplot(222)
    p2=plot(q_mt*180/pi,h*1e3,'--','DisplayName',[num2str(order) '^e order polynomial']);
    
    rel_err2 = abs((h-MA_PF)./MA_PF);
    err2 = max(rel_err2);
    subplot(224)
    semilogy(q_mt*180/pi,rel_err2,'Color',p2.Color)
    hold on

%     disp(err2)
    
    if err2 <= 1e-5 || order >= 5
        stop = 1;
    else
        order = order+1;
    end
end

subplot(224)
hold on
grid on
xlabel('Midtarsal angle (°)')
ylabel('log_{10}\deltaMA_{PF} (-)')
title('Relative error of the approximation')
axis tight

codeString2 = ['MA_PF = ' num2str(coeff2(1),7)];
for ii=1:length(coeff2)-1
    codeString2 = [codeString2 ' + ' num2str(coeff2(ii+1),7) '.*q_mt.^' num2str(ii)];
end
codeString2 = [codeString2 ';'];
disp(codeString2)


%% implementing the generated code, to see there are enough decimals

l_PF_fa_p = 0.1392179 + 0.0374482.*q_mt.^1 + -0.0166876.*q_mt.^2 + -0.001758651.*q_mt.^3 + 0.0004480769.*q_mt.^4;
MA_PF_p = 0.0374478 + -0.03337403.*q_mt.^1 + -0.005255987.*q_mt.^2 + 0.001767266.*q_mt.^3 + -0.0001071423.*q_mt.^4 + 9.858065e-05.*q_mt.^5;

subplot(221)
p1=plot(q_mt*180/pi,l_PF_fa_p*1e3,'DisplayName','code');

subplot(223)
semilogy(q_mt*180/pi,abs((l_PF_fa_p-l_PF_fa)./l_PF_fa),'Color',p1.Color)

subplot(222)
p2=plot(q_mt*180/pi,MA_PF_p*1e3,'DisplayName','code');

subplot(224)
semilogy(q_mt*180/pi,abs((MA_PF_p-MA_PF)./MA_PF),'Color',p2.Color)

%%

% figure
% plot(q_mt*180/pi,MA_PF*1e3,'DisplayName','exact')
% hold on
% grid on
% % legend('Location','northeast');
% xlabel('Midtarsal angle (°)')
% ylabel('MA_{PF,fa} (mm)')
% title('Moment arm PF')


%%
% order = 5;
% 
% [mat,diff_mat_q] = n_art_mat_3(q_mt, order);
% 
% coeff=[mat ; diff_mat_q]\[l_PF_fa; MA_PF];
% dM_recon = diff_mat_q*coeff;
% lMT_recon=mat*coeff;
% 
% subplot(221)
% plot(q_mt*180/pi,lMT_recon*1e3)
% 
% subplot(222)
% plot(q_mt*180/pi,dM_recon*1e3)
