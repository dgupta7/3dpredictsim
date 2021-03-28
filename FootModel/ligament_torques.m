bounds_mtj = [-35*0.75;35]*pi/180;


q_mt = linspace(-35,35,500)'*pi/180;

k_pass.mtj = [-30 5 5 -3]';
theta.pass.mtj = bounds_mtj;
tau_pass = k_pass.mtj(1,1)*exp(k_pass.mtj(2,1)*(q_mt-theta.pass.mtj(2,1))) + ...
    k_pass.mtj(3,1)*exp(k_pass.mtj(4,1)*(q_mt-theta.pass.mtj(1,1)));

q0 = 5*pi/180;
kMT_li = 10;
M_li = -(q_mt - q0)*kMT_li;
M_li = M_li + tau_pass;

figure
plot(q_mt,M_li)
grid on
