function main2()
model = getModel(2);

endTime = 1000*10^-3;
timeStep = 10^-3;
i=1;
for t=0:timeStep:endTime
    [qdd,r] = LE_solver(model);
    R1(i,:) = r(1:3,1)';
    R2(i,:) = r(1:3,2)';
    [q,qd] = motionIntegration(model.q,model.qd,qdd,timeStep);
    %updating model
    model.q = q;
    model.qd = qd;
    %simulateChain(model);
    i=i+1;
end

M = dlmread('twoLink3D_cm_pos.tab','\t');
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
title('link-1x')
legend('lag','adams')
figure()
hold on
grid on
plot(t,R1(:,2),'g',t,r1(:,2),'k');
title('link-1y')
legend('lag','adams')
figure()
hold on
grid on
plot(t,R1(:,3),'b',t,r1(:,3),'m');
title('link-1z')
legend('lag','adams')

figure()
hold on
grid on
plot(t,R2(:,1),'r',t,r2(:,1),'c');
title('link-2x')
legend('lag','adams')
figure()
hold on
grid on
plot(t,R2(:,2),'g',t,r2(:,2),'k');
title('link-2y')
legend('lag','adams')
figure()
hold on
grid on
plot(t,R2(:,3),'b',t,r2(:,3),'m');
title('link-2z')
legend('lag','adams')


end