

n = 500;

q_mtj = linspace(-15,15,n)*pi/180;


M_li = 0.081222 + -10.873*q_mtj.^1 + -307.95*q_mtj.^2 + -3487.8*q_mtj.^3 + -7202.7*q_mtj.^4 + 76515*q_mtj.^5;
figure
plot(q_mtj*180/pi,M_li)
hold on

%%

K_pass = [-60 30 60 -30]';
theta_pass = [-0.3 0.3]';

Tau_pass = K_pass(1,1)*exp(K_pass(2,1)*(q_mtj-theta_pass(2,1))) + ...
    K_pass(3,1)*exp(K_pass(4,1)*(q_mtj-theta_pass(1,1)));

plot(q_mtj*180/pi,Tau_pass,'--')


