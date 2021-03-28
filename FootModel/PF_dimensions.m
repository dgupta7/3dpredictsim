% plantar fascia dimensions

% calcaneal tubercle
w = 15.45;
t = 2.79;
l_1 = 146.45;
c1 = w*t/l_1;

% central aponeurosis bundle
w1 = 71.2;
t1 = 1.26;
l1 = 42.12;
w2 = 6.68;
t2 = 1.04;
l2 = 44.43;
w3 = 6.17;
t3 = 0.91;
l3 = 41.43;
w4 = 4.93;
t4 = 0.84;
l4 = 36.48;
w5 = 4.55;
t5 = 0.72;
l5 = 40.09;

c2 = w1*t1/l1 + w2*t2/l2 + w3*t3/l3 + w4*t4/l4 + w5*t5/l5;

A0 = c1*c2/(c1+c2) *(l_1+l3);

