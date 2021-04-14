% When keeping the toes at a constant angle relative to the ground, changes
% in midtarsal angle affect the mtp angle.

calcn2mtj = 0.08207;
mtj2mtpj = 0.089638;
beta0 = 2.4935;

q_mt = 0*pi/180;

beta = beta0 + q_mt; % top angle of foot arch triangle
l_fa = sqrt(calcn2mtj^2 + mtj2mtpj^2 - 2*calcn2mtj*mtj2mtpj*cos(beta)); % foot arch length
h_fa = calcn2mtj*mtj2mtpj/l_fa*sin(beta); % foot arch height

gamma0 = asin(h_fa/mtj2mtpj);


qs_mt = linspace(-35,35,20)';
qs_mtp = [-45:15:45]';

for j=1:length(qs_mtp)
    for i=1:length(qs_mt)

        q_mtp0 = qs_mtp(j);
        q_mt = qs_mt(i)*pi/180;

        beta = beta0 + q_mt; % top angle of foot arch triangle
        l_fa = sqrt(calcn2mtj^2 + mtj2mtpj^2 - 2*calcn2mtj*mtj2mtpj*cos(beta)); % foot arch length
        h_fa = calcn2mtj*mtj2mtpj/l_fa*sin(beta); % foot arch height

        gamma = asin(h_fa/mtj2mtpj);

        dq(i,j) = (gamma-gamma0)*180/pi;
        q_mtp(i,j) = q_mtp0 + dq(i,j);
        
    end
end

figure
plot(qs_mt,q_mtp)

figure
plot(qs_mt,dq)

c = nanmean(dq./qs_mt);

q_mtp = qs_mtp' + c.*qs_mt;
