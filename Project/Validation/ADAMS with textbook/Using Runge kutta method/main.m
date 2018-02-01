function main(NB)
model = getModel(NB);

endTime = 5000*10^-3;
timeStep = 10^-3;
i=1;
for t=0:timeStep:endTime
    [qdd,r] = EX2_solver(model);
    R1(i,:) = r(1:3,1)';
    R2(i,:) = r(1:3,2)';
    
    q = model.q;
    qd = model.qd;
    q1 = q;
    qd1 = qd;
    qdd1 = qdd;
    
    q2 = q + 0.5*qd1*timeStep;
    qd2 = qd + 0.5*qdd1*timeStep;
    model.q = q2;
    model.qd = qd2;
    qdd2 = EX2_solver(model);

    q3 = q + 0.5*qd2*timeStep;
    qd3 = qd + 0.5*qdd2*timeStep;
    model.q = q3;
    model.qd = qd3;
    qdd3 = EX2_solver(model);

    q4 = q + qd3*timeStep;
    qd4 = qd + qdd3*timeStep;
    model.q = q4;
    model.qd = qd4;
    qdd4 = EX2_solver(model);

    qf = q + (timeStep/6.0)*(qd1 + 2*qd2 + 2*qd3 + qd4);
    qdf = qd + (timeStep/6.0)*(qdd1 + 2*qdd2 + 2*qdd3 + qdd4);
    
    %updating model
    model.q = qf;
    model.qd = qdf;
    %simulateChain(model);
    i=i+1;
end
M = dlmread('twoExample_5.tab','\t');
r1 = M(:,2:4)./1000;
r2 = M(:,5:7)./1000;
t=0:timeStep:endTime;

%calculaing maximum error
diff1 = R1-r1;
diff2 = R2-r2;
[C,I] = max(abs(diff1));
[C,I1] = max(C);
C,I(I1)
[C,I] = max(abs(diff2));
[C,I1] = max(C);
C,I(I1)

figure()
hold on
grid on
plot(t,R1(:,1),'r',t,r1(:,1),'c');
plot(t,R1(:,2),'g',t,r1(:,2),'k');
plot(t,R1(:,3),'b',t,r1(:,3),'m');
title('Link 1')
legend('x-lag','x-adams','y-lag','y-adams','z-lag','z-adams')

figure()
hold on
grid on
plot(t,R2(:,1),'r',t,r2(:,1),'c');
plot(t,R2(:,2),'g',t,r2(:,2),'k');
plot(t,R2(:,3),'b',t,r2(:,3),'m');
title('Link 2')
legend('x-lag','x-adams','y-lag','y-adams','z-lag','z-adams')
end

