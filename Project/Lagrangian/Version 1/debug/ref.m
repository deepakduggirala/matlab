g=9.806;
tht = 0;
omega = 0;
l=0.1;
delT = 10^-3;
figure()
axis equal
hold on
grid on
xlabel('X-axis')
ylabel('Y-axis')
for t=0:delT:0.2
    alpha = 3*g*cos(tht)/(2*l);
    omega = omega + alpha*delT;
    tht = tht + omega*delT;
    x =l*cos(tht);
    y = -l*sin(tht);
    plot(x,y,'o')
end

