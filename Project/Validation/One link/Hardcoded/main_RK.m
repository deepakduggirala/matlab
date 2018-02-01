function main_RK()
q=0;
qd=0;
l = 0.1;
endTime = 5000*10^-3;
timeStep = 10^-3;
i=1;
g = -9.81;
for t=0:timeStep:endTime
    qdd = 1.5*g*cos(q)/l;
    R(i)=qdd;
    
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

    
    i=i+1;
end
M = dlmread('oneLink_5_0.001_less_error_ang_acc.tab','\t');
r1 = degtorad(M(:,2));

t=0:timeStep:endTime;
%calculaing maximum error
diff1 = R'-r1;
% diff2 = R2-r2;
[C,I] = max(abs(diff1));
[C,I1] = max(C);
C,I(I1)

grid on
plot(t,R,'g',t,r1,'k');
title('Link 1')
legend('ang_acc-text','ang_acc-adams')
end