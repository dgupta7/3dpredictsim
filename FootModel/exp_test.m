

c1 = 1;
c2 = 10;
c3 = 0.1;

syms 'x'

y(x) = c1*exp(c2*(x-c3));

dy(x) = diff(y,x);


xs = linspace(-10,20,500)'*pi/180;

for i=1:length(xs)
    ys(i) = y(xs(i));
    
end

k0 = dy(c3);

xt = [c3-0.1;c3;c3+0.1];
yt = c1+[-0.1;0;0.1]*k0;


figure
plot(xs*180/pi,ys)
hold on
plot(xt*180/pi,yt)

