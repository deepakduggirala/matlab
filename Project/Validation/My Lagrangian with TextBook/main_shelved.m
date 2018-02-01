function main(NB)
model_TB = getModel(NB);
model_LG = getModel(NB);

endTime = 5000*10^-3;
timeStep = 10^-3;
i=1;
fileId_TB = fopen('TB.txt','w');
fileId_LG = fopen('LG.txt','w');

for t=0:timeStep:endTime
    [qdd_TB,r_TB] = EX2_solver(model_TB);
    [qdd_LG,r_LG] = LE_solver(model_LG);
    
    printToFile2(fileId_TB,t,r_TB);
    printToFile2(fileId_LG,t,r_LG);
    
    R1_TB(i,:) = r_TB(1:3,1)';
    R1_LG(i,:) = r_LG(1:3,1)';
    R2_TB(i,:) = r_TB(1:3,2)';
    R2_LG(i,:) = r_LG(1:3,2)';
    [q_TB,qd_TB] = motionIntegration(model_TB.q, model_TB.qd ,qdd_TB, timeStep);
    [q_LG,qd_LG] = motionIntegration(model_LG.q, model_LG.qd ,qdd_LG, timeStep);
    %updating model
    model_TB.q = q_TB;
    model_LG.q = q_LG;
    model.qd_TB = qd_TB;
    model.qd_LG = qd_LG;
    %simulateChain(model);
    
    i=i+1;
end

%calculaing maximum error
diff1 = R1_TB-R1_LG;
diff2 = R2_TB-R2_LG;
[C,I] = max(abs(diff1));
[C,I1] = max(C);
C,I(I1)
[C,I] = max(abs(diff2));
[C,I1] = max(C);
C,I(I1)

t=0:timeStep:endTime;

figure()
hold on
grid on
plot(t,R1_TB(:,1),'r',t,R1_LG(:,1),'c');
plot(t,R1_TB(:,2),'g',t,R1_LG(:,2),'k');
plot(t,R1_TB(:,3),'b',t,R1_LG(:,3),'m');
title('Link 1')
legend('x-TB','x-LG','y-TB','y-LG','z-TB','z-LG')

figure()
hold on
grid on
plot(t,R2_TB(:,1),'r',t,R2_LG(:,1),'c');
plot(t,R2_TB(:,2),'g',t,R2_LG(:,2),'k');
plot(t,R2_TB(:,3),'b',t,R2_LG(:,3),'m');
title('Link 2')
legend('x-TB','x-LG','y-TB','y-LG','z-TB','z-LG')
end

