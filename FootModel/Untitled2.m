
% model 0 dof
m0.t.mtpj = m0.t.subt + m0.cmf.mtpj;
m0.t.cmfCOM = m0.t.subt + m0.cmf.COM;

m0.t.cmfI = Steiner(-m0.t.cmfCOM,m0.cmf.I,m0.cmf.m);

m0.c.subt_st = [0.78717961, 0.60474746, -0.12094949]';
m0.t.subt_st = m0.t.subt + m0.c.subt_st/10;

m0j = [[0;0;0],m0.t.subt,m0.t.mtpj];
m0c = m0.t.cmfCOM;
m0st = [m0.t.subt,m0.t.subt_st];


% model 2 dof
m2.t.mtj = m2.t.subt + m2.c.mtj;
m2.t.tmtj = m2.t.mtj + m2.m.tmtj;
m2.t.mtpj = m2.t.tmtj + m2.f.mtpj;

m2.c.subt_st = [-0.63276, -0.588866, 0.502838]';
m2.t.subt_st = m2.t.subt + m2.c.subt_st/10;

m2.t.cCOM = m2.t.subt + m2.c.COM;
m2.t.mCOM = m2.t.mtj + m2.m.COM;
m2.t.fCOM = m2.t.tmtj + m2.f.COM;

m2.t.cmfCOM = (m2.t.cCOM*m2.c.m + m2.t.mCOM*m2.m.m + m2.t.fCOM*m2.f.m)/(m2.c.m+m2.m.m+m2.f.m);

m_sf = m0.cmf.m/(m2.c.m+m2.m.m+m2.f.m);

m2j = [[0;0;0],m2.t.subt,m2.t.mtj,m2.t.tmtj,m2.t.mtpj];
m2c = [m2.t.cCOM,m2.t.mCOM,m2.t.fCOM];
m2st = [m2.t.subt,m2.t.subt_st];

sf = m0.gnd.talus2mtpj./m2.gnd.talus2mtpj;
sf(2) = 1;

m2js(1,:) = m2j(1,:)*sf(1);
m2js(2,:) = m2j(2,:)*sf(2);
m2js(3,:) = m2j(3,:)*sf(3);
m2cs(1,:) = m2c(1,:)*sf(1);
m2cs(2,:) = m2c(2,:)*sf(2);
m2cs(3,:) = m2c(3,:)*sf(3);
m2.t.cmfCOM_s = m2.t.cmfCOM.*sf;

% model 1 with mtj:
diff = m0.t.mtpj(2) - m2.t.mtpj(2);
m11.t.subt = m0.t.subt;
m11.t.subt(2) = m11.t.subt(2) - diff;
m11.t.mtj = m2.t.mtj.*sf;
m11.t.mtpj = m0.t.mtpj;
m11.t.mtpj(2) = m11.t.mtpj(2) - diff;

m11.c.mtj = m11.t.mtj - m11.t.subt;
m11.mf.mtpj = m11.t.mtpj - m11.t.mtj;

m11.t.cCOM = m2.t.cCOM.*sf;
m11.t.cCOM(2) = m11.t.cCOM(2) + 0.01;
m11.t.mfCOM = (m2.t.mCOM*m2.m.m + m2.t.fCOM*m2.f.m)/(m2.m.m+m2.f.m).*sf;
m11.t.cmfCOM = m0.t.cmfCOM;
m11.t.cmfCOM(2) = m2.t.cmfCOM_s(2);

m11.c.m = m2.c.m*m_sf;
m11.mf.m = (m2.m.m+m2.f.m)*m_sf;

m11.c.COM = m11.t.cCOM - m11.t.subt;
m11.mf.COM = m11.t.mfCOM - m11.t.mtj;

m11j = [[0;0;0],m11.t.subt,m11.t.mtj,m11.t.mtpj];
m11c = [m11.t.cCOM,m11.t.mfCOM];

Icx = 0.08842e-3;     %kg m^2 /kg (normalised with cmf mass)
Icy = 0.1807e-3;
Icz = 0.2077e-3;
m11.c.I = [Icx; Icy; Icz]*m0.cmf.m;

m11.t.cI = Steiner(-m11.t.cCOM,m11.c.I,m11.c.m);

m11.mf.I = [0;0;0]; %initialise
m11.t.mfI = Steiner(-m11.t.mfCOM,m11.mf.I,m11.mf.m);


m0.t.cmfI = Steiner(-m11.t.cmfCOM,m0.cmf.I,m0.cmf.m); % Steiner

I_res = m0.t.cmfI - (m11.t.cI + m11.t.mfI); % resultant error is I_mf value to satisfy
m11.mf.I = I_res;

locSphere_3_r_new = locSphere_3_r - m11.c.mtj';
locSphere_4_r_new = locSphere_4_r - m11.c.mtj';

m11.mf.ls3 = locSphere_3_r_new';
m11.mf.ls4 = locSphere_4_r_new';

m11.t.ls1 = m11.t.subt + locSphere_1_r';
m11.t.ls2 = m11.t.subt + locSphere_2_r';
m11.t.ls3 = m11.t.mtj + m11.mf.ls3;
m11.t.ls4 = m11.t.mtj + m11.mf.ls4;
m11.t.ls5 = m11.t.mtpj + locSphere_5_r';
m11.t.ls6 = m11.t.mtpj + locSphere_6_r';

locSps = [m11.t.ls1,m11.t.ls2,m11.t.ls3,m11.t.ls4,m11.t.ls5,m11.t.ls6];

figure
subplot(2,1,1)
hold on
grid on

plot(m2j(1,:),m2j(2,:),'ob')
plot(m2c(1,:),m2c(2,:),'xb')
plot(m2.t.cmfCOM(1),m2.t.cmfCOM(2),'*b')
plot(m2st(1,:),m2st(2,:),'--b')

plot(m2js(1,:),m2js(2,:),'ok')
plot(m2js(1,:),m2js(2,:),'.k')
plot(m2cs(1,:),m2cs(2,:),'xk')
plot(m2.t.cmfCOM_s(1),m2.t.cmfCOM_s(2),'*k')

plot(m0j(1,:),m0j(2,:),'og')
plot(m0c(1,:),m0c(2,:),'xg')
plot(m0st(1,:),m0st(2,:),'--g')

plot(m11j(1,:),m11j(2,:),'or')
plot(m11c(1,:),m11c(2,:),'xr')

plot(locSps(1,:),locSps(2,:),'*c')
axis equal

subplot(2,1,2)
hold on
grid on

plot(m2j(1,:),m2j(3,:),'ob')
plot(m2c(1,:),m2c(3,:),'xb')
plot(m2.t.cmfCOM(1),m2.t.cmfCOM(3),'*b')
plot(m2st(1,:),m2st(3,:),'--b')

plot(m2js(1,:),m2js(3,:),'ok')
plot(m2js(1,:),m2js(3,:),'.k')
plot(m2cs(1,:),m2cs(3,:),'xk')
plot(m2.t.cmfCOM_s(1),m2.t.cmfCOM_s(3),'*k')

plot(m0j(1,:),m0j(3,:),'og')
plot(m0c(1,:),m0c(3,:),'xg')
plot(m0st(1,:),m0st(3,:),'--g')

plot(m11j(1,:),m11j(3,:),'or')
plot(m11c(1,:),m11c(3,:),'xr')

plot(locSps(1,:),locSps(3,:),'*c')
axis equal

%% relating vectors to foot arch compression
% based on DOI: 10.1038/srep19403

a = m11.c.mtj(1:2);
b = m11.mf.mtpj(1:2);

l_0 = norm(a+b);
h0 = -b(2); 
% Difference in y-coordinate between calcn origin and m1a.f.mtpj 
% results in 1° difference, so it is omitted.

c0 = acos(h0/norm(a));
d0 = acos(h0/norm(b));

mt0 = (c0+d0)*180/pi;

h1 = h0*0.8; 
c1 = acos(h1/norm(a));
d1 = acos(h1/norm(b));

mt1 = (c1+d1)*180/pi;

mt_bound = mt1 - mt0;
% So 15° bound is sensible

h2 = h0*0.87; 
c2 = acos(h2/norm(a));
d2 = acos(h2/norm(b));

mt2 = (c2+d2)*180/pi;

mt_2 = mt2 - mt0;

%% subtalar joint axis orientation
% original
alpha_st = atan(m0.c.subt_st(2)/m0.c.subt_st(1))*180/pi;
beta_st = atan(-m0.c.subt_st(3)/m0.c.subt_st(1))*180/pi;

incl1 = alpha_st*pi/180;
dev1 = beta_st*pi/180;

R1 = [cos(incl1),-sin(incl1),0;
     sin(incl1),cos(incl1),0;
     0,0,1];
R2 = [cos(dev1),0,sin(dev1);
     0,1,0;
     -sin(dev1),0,cos(dev1)];

st_subt1 = R2*R1*[1;0;0]

m0.c.subt_st;

% doi:10.1136/bjsm.2010.080119
incl2 = 42*pi/180; % +-16
dev2 = 11*pi/180; % +-23

R1 = [cos(incl2),-sin(incl2),0;
     sin(incl2),cos(incl2),0;
     0,0,1];
R2 = [cos(dev2),0,sin(dev2);
     0,1,0;
     -sin(dev2),0,cos(dev2)];

st_subt2 = R2*R1*[1;0;0]
norm(st_subt2);

% doi:10.1016/j.jbiomech.2012.01.011
incl3 = 45.5*pi/180; % std 3.4
dev3 = 5*pi/180; % std3.8


R1 = [cos(incl3),-sin(incl3),0;
     sin(incl3),cos(incl3),0;
     0,0,1];
R2 = [cos(dev3),0,sin(dev3);
     0,1,0;
     -sin(dev3),0,cos(dev3)];

st_subt3 = R2*R1*[1;0;0]
norm(st_subt3);






%%
function I_new = Steiner(vec1,I_com,mass)
    I_new = I_com + mass*(vec1'*vec1*eye(3) - vec1*vec1')*ones(3,1);

end

