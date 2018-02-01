function proj_v0_0_42_debug()
global r0;
global r;
global J;
global J0;
global m;

p=1;
elemLen = 0.1;
m = p*elemLen;
endTime = 0.2;
timeStep = 10^-3;

tau = 0;
thtd=0;
tht=degtorad(0);

r0 = zeros(3,1);
r0(1,:) = -0.5*elemLen;
r = r0;

J0 = [m*(elemLen^2)/12 0 0;
    0 0 0;
    0 0 0];
J=J0;
figure()
axis equal
grid on
hold on
xlabel('X-axis')
ylabel('Y-axis')
for t=0:timeStep:endTime
    myPlot(tht)
    D = getD(tht);
    H = 0;
    C = getC(tht);

    thtdd = (tau-H-C)/D;
    thtd = thtd + thtdd*timeStep;
    tht = tht + thtd*timeStep;
    update_r_J(tht);
end
end

function D= getD(tht)
global J
U11 = [-sin(tht) cos(tht) 0;
       -cos(tht) -sin(tht) 0;
       0 0 0];
   D = trace(U11*J*(U11'));
end

function C = getC(tht)
global r0
global m
U11 = [-sin(tht) cos(tht) 0;
       -cos(tht) -sin(tht) 0;
       0 0 0];
C = -m*[0 -9.8 0]*U11*r0;
end

function update_r_J(tht)
global r0;
global r;
global J;
global J0;
    R = [cos(tht) sin(tht) 0;
        -sin(tht) cos(tht) 0;
        0 0 1];
    r = R*r0;
    J =R*J0*(R');
end

function myPlot(tht)
    y = 0.1*sin(tht);
    x = 0.1*cos(tht);
    plot(x,y,'o')
end
