function main3(NB)
model = getModel(NB);

endTime = 10*10^-3;
timeStep = 10^-3;
i=1;
for t=0:timeStep:endTime
    [qdd,r,H,H1] = LE_solver(model);
    H-H1
    diff(:,i) = H-H1;
    [q,qd] = motionIntegration(model.q,model.qd,qdd,timeStep);
    %updating model
    model.q = q;
    model.qd = qd;
    %simulateChain(model);
    i=i+1;
end
[C,I] = max(abs(diff));
[C,I1] = max(C);
C,I(I1)
end