

% vars_init = [-10*pi/180;0.45;10*pi/180;0;0]./scale_qs;
% [res] = f_foot([vars_init;zeros(NMf,1)],0,0);

% full(res)
% full(res(18:end))

[res2,tm1,tm2,tm3,tm4] = f_foot_pp([qs_sol./scale_qs;FTs_sol./scale_FTs],0,0);
full(res2)
full(tm1)
full(tm2)
full(tm3)
full(tm4)


%%

[Test] = full(F(zeros(30,1)));
toes_or = Test(jointfi.toes_or);
metatarsi_or = Test(jointfi.metatarsi_or);
calcn_or = Test(jointfi.calcn_or);
talus_or = Test(jointfi.talus_or);
tibia_or = Test(jointfi.tibia_or);

x = [toes_or(1), metatarsi_or(1), calcn_or(1), talus_or(1), tibia_or(1)];
y = [toes_or(2), metatarsi_or(2), calcn_or(2), talus_or(2), tibia_or(2)];
z = [toes_or(3), metatarsi_or(3), calcn_or(3), talus_or(3), tibia_or(3)];
figure
subplot(2,1,1)
plot(x,y,'-o')
axis equal
subplot(2,1,2)
plot(x,z,'-o')
axis equal

    
% l_fa_ext = norm(squeeze(toes_or(1:2)-calcn_or(1:2)));
% % h_fa_ext = metatarsi_or(2)-calcn_or(2);
% x0 = metatarsi_or(1);
% x1 = calcn_or(1);
% x2 = toes_or(1);
% y0 = metatarsi_or(2);
% y1 = calcn_or(2);
% y2 = toes_or(2);
% 
% h = abs( (x2-x1)*(y1-y0) - (x1-x0)*(y2-y1) )/(sqrt( (x2-x1)^2 + (y2-y1)^2 ) )



