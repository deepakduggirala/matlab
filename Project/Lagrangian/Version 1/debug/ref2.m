function ref2()
m1 = 1;
m2 = 1;

l = 1;
tht1=0;
tht2=0;
thtd1=0;
thtd2=0;
g=9.8;
nL = 2;
endTime = 10*10^-3;
timeStep = 10^-3;
thtd=zeros(nL,1);
tht=zeros(nL,1);
for t=0:timeStep:endTime
    D = (10^-3)*[(5/3)+cos(tht2) (1/3)+cos(tht2)/2; (1/3)+cos(tht2)/2  (1/3)]
    H = (10^-3)*[-0.5*sin(tht2)*(thtd2)^2 - sin(tht2)*thtd1*thtd2; 0.5*sin(tht2)*(thtd1^2)]
    C = (10^-2)*9.806*[0.5*cos(tht1) + 0.5*cos(tht1+tht2) + cos(tht1); 0.5*cos(tht1 + tht2)]
    
    thtdd = inv(D)*([0;0]-H-C)
    thtd = thtd + thtdd*timeStep
    tht = tht + thtd*timeStep
end
end
