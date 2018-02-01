function main()
q=0;
qd=0;
l = 0.1;
endTime = 5000*10^-3;
timeStep = 10^-3;
i=1;
g = -9.81;
for t=0:timeStep:endTime
    qdd = 1.5*g*cos(q)/l;
    qd = qd+qdd*timeStep;
    q = q+qd*timeStep;
    R(i)=qdd;
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