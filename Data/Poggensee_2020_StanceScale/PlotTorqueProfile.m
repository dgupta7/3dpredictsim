%% Plot torque profile

load('torque_profile.mat');
lw = 2;

figure();
plot(time,torque,'k','LineWidth',lw);
xlabel('% gait cycle');
ylabel('Ankle torque [Nm]');
delete_box