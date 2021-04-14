%% old models for the midtarsal torque resulting from passive structures excluding plantar fascia

% dl_0 = 3*pi/180;
% dq_0 = 2*pi/180;
% M_li = -kMT_li*( (q_mt-dq_0) - dl_0*tanh((q_mt-dq_0)/dl_0/1.2)) +... % toe-in and linear part
%         (-exp(20*(q_mt-20*pi/180)) + exp(-25*(q_mt+10*pi/180)))/2; % stiffening at the end
% 
% M_li = -9*(exp(4*(q_mt-2*pi/180))-1)*1.2 + 2*exp(-10*(q_mt+0.1));
% 
% M_li = -kMT_li*q_mt*(exp(5*(q_mt-10*pi/180)) + exp(-3*(q_mt+20*pi/180)))/2; % approx 1
    

% k1 = 12;
% t1 = 8*pi/180*sf_li;
% f2 = 1.5*(2-sf_li);
% k2 = 5;
% t2 = 10*pi/180*sf_li;
% M_li = (-exp(k1*(q_mt-t1)) + f2*exp(-k2*(q_mt+t2)))*kMT_li*5/90; % approx 2
% M_li = (-exp(k1*(q_mt-t1)) + f2*exp(-k2*(q_mt+t2)))*5*3-2; % approx 2a
% 
% sf = 0.1;
% t = 0.5;
% c1 = 12*(2-t)/(1-t+sf);
% c2 = 8*sf;
% c3 = 1.5;
% c4 = 5;
% c5 = 10*sf;
% M_li = (-exp(c1*(q_mt-c2*pi/180)) + c3*exp(-c4*(q_mt+c5*pi/180)))*kMT_li*5/90*(1+sf)/2;
% 
% sf = 0.3;
% c1 = 12/(sf+0.5)*1.5;
% c2 = 8*(1+sf)/2;
% c3 = 1.5*(sf);%+0.4)/1.4;
% c4 = 5;%/(sf+0.5)*1.5;
% c5 = 10;%*(1+sf)/2;

% sf = 0.1;
% c1 = 12/(sf+1)*2;
% c2 = 8*(2+sf)/3;
% c3 = 1.5*(sf+0.5)/1.5;
% c4 = 5;%/(sf+5)*6;
% c5 = 10;%*(1+sf)/2;

% c1 = 25;
% c2 = 1;
% c3 = 3;
% c4 = 5;
% c5 = 1;
% sf=1;
% 
% M_li = (-exp(c1*(q_mt-c2*pi/180)) + c3*exp(-c4*(q_mt+c5*pi/180)))*5/sf;

% M_li = (-exp(12*(q_mt*3-8*pi/180)) + 2*exp(-5*(q_mt*2+10*pi/180)))*5*2;

% M_li = (-exp(12*(q_mt*2-8*pi/180)) + 1.5*exp(-5*(q_mt*1.5+10*pi/180)))*5*5; % v5

% M_li = 0;

% M_li = -exp(25*(q_mt-5*pi/180)) + 2*exp(-15*(q_mt+5*pi/180)) +0.5;

% M_li = -2*exp(10*(q_mt-5*pi/180)) + 2*exp(-15*(q_mt+5*pi/180));

% M_li = -19*exp(3*(q_mt-10*pi/180)) + 19*exp(-3*(q_mt+10*pi/180)); % Ker

