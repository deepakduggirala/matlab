function mainQ(NB)
timeStep = 10^-3;
endTime  = 1000*10^-3;

t_ref = 10^-3;
e_ref = 1000*10^-3;

refDataFileName = 'L3_3D_Q.tab';

model = getModel1(NB);

i=1;
for t=0:timeStep:endTime
    qdd = LE_solver(model);
    model.qd = model.qd + qdd*timeStep;
    model.q = model.q + model.qd*timeStep;
    Q(i,:) = rad2deg(model.q');
    i=i+1;
end
t1 = 0:timeStep:endTime;
t2 = 0:t_ref:e_ref;

if(strcmp(refDataFileName,'') == 0)
    flag = 1;
    M = dlmread(refDataFileName,'\t');
end

H1 = figure();
grid on;
hold on;
plot(t1,Q(:,1),'g');
if(flag)
    plot(t2,M(:,2),'r');
end
title('link1');
legend('LAG','ADAMS');

H2 = figure();
grid on;
hold on;
plot(t1,Q(:,2),'g');
if(flag)
    plot(t2,M(:,3),'r');
end
title('link2');
legend('LAG','ADAMS')


H3 = figure();
grid on;
hold on;
plot(t1,Q(:,3),'g');
if(flag)
    plot(t2,M(:,4),'r');
end
title('link3');
legend('LAG','ADAMS')
end