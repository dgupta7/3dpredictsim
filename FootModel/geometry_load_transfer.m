

a = 0.08207;
b = 0.089638;
phi0 = 2.493499;
H0 = 0.027280;

l1 = sqrt(a^2-H0^2);
l2 = sqrt(b^2-H0^2);

l = l1 + l2;

c = 0.03345;
h = H0;

h = linspace(0.02,0.05,100)';

fr1 = ( (l1*(l2+c))/l - c )./h; %F/load for back side
fr2 = ( l2*(l1-c)/l )./h; %F/load for front side

figure
plot(h*1e3,fr1)
xlabel('resting arch height (mm)')
ylabel('PF force / vert load on talus (-)')
title('Force trasfer talus to PF in initial position')


%%
h = H0;

c = linspace(0,0.05,100)';

fr1 = ( (l1*(l2+c))/l - c )./h; %F/load for back side

figure
plot(c*1e3,fr1)
xlabel('talus to mtj distance (mm)')
ylabel('PF force / vert load on talus (-)')
title('Force trasfer talus to PF in initial position')


%%


