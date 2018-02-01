function main(NB)
    model = getModel(NB);
    endTime = 1000*10^-3;
    timeStep = 10^-3;
    
    i=1;
    for t=0:timeStep:endTime
        qdd = LE_solver(model);
        model.qd = model.qd + qdd*timeStep;
        model.q = model.q + model.qd*timeStep;
        Q(i,:) = model.q';
        i=i+1;
    end
    t1 = 0:timeStep:endTime;
    hX = figure();
    hold on
    grid on
    plot(t1,Q(:,1),'r');
    title('1');
    
    hY = figure();
    hold on
    grid on
    plot(t1,Q(:,2),'r')
    title('2');
    
    hZ = figure();
    hold on
    grid on
    plot(t1,Q(:,3),'r')
    title('3');
end