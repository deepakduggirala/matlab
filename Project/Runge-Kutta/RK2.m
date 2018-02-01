function RK2
endTime = 1000*10^-3;
delT = 10^-3;
timeStep = delT;
m_RK = getModel(2);
m_EU = getModel(2);

t_ref = 10^-3;
e_ref = 1000*10^-3;

refDataFileName = 'twoLink2DQ.tab';

i=1;
for t=0:delT:endTime
    m    = m_RK;
    q1   = m.q;
    qd1  = m.qd;
    qdd1 = LE_solver(m);
    
    q2   = q1  + 0.5*qd1*delT + qdd1*(delT^2)/8;
    qd2  = qd1 + 0.5*qdd1*delT;
    m.q  = q2;
    m.qd = qd2;
    qdd2 = LE_solver(m);
    
    q3   = q1  + 0.5*qd2*delT + qdd2*(delT^2)/8;
    qd3  = qd1 + 0.5*qdd2*delT;
    m.q  = q3;
    m.qd = qd3;
    qdd3 = LE_solver(m);
    
    q4   = q1  + qd3*delT + qdd3*(delT^2)/2;
    qd4  = qd1 + qdd3*delT;
    m.q  = q4;
    m.qd = qd4;
    qdd4 = LE_solver(m);
    
    qf   = q1  + (delT/6)*(qd1 + 2*qd2 + 2*qd3 + qd4);
    qdf  = qd1 + (delT/6)*(qdd1 + 2*qdd2 + 2*qdd3 + qdd4);
    m_RK.q = qf;
    m_RK.qd = qdf;
    
    Q_RK(i,:) = rad2deg(m_RK.q');
    i=i+1;
end

i=1;
for t = 0:delT:endTime
    qdd = LE_solver(m_EU);
    qd_start = m_EU.qd;
    m_EU.qd = m_EU.qd + qdd*delT;
    m_EU.q = m_EU.q + qd_start*delT + qdd*(delT^2)/2;
    
    Q_EU(i,:) = rad2deg(m_EU.q');
    i=i+1;
end

t1 = 0:timeStep:endTime;
t2 = 0:t_ref:e_ref;

if(strcmp(refDataFileName,'') == 0)
    flag = 1;
    M = dlmread(refDataFileName,'\t');
end
%calculation of maximum absolute deviation
max(abs(Q_EU(:,1) - M(:,2)))
max(abs(Q_EU(:,2) - M(:,3)))
max(abs(Q_RK(:,1) - M(:,2)))
max(abs(Q_RK(:,2) - M(:,3)))



H1 = figure();
grid on;
hold on;
plot(t1,Q_EU(:,1),'g');
if(flag)
    plot(t2,M(:,2),'r');
end
plot(t1,Q_RK(:,1),'b');
title('link1');
legend('EU','ADAMS','RK');

H2 = figure();
grid on;
hold on;
plot(t1,Q_EU(:,2),'g');
if(flag)
    plot(t2,M(:,3),'r');
end
plot(t1,Q_RK(:,2),'b');
title('link2');
legend('EU','ADAMS','RK');


% H3 = figure();
% grid on;
% hold on;
% plot(t1,Q_EU(:,3),'g');
% if(flag)
%     plot(t2,M(:,4),'r');
% end
% plot(t1,Q_RK(:,3),'b');
% title('link3');
% legend('EU','ADAMS','RK');

H4 = figure();
grid on
hold on
plot(t1,Q_EU(:,1)-M(:,2),'g')
plot(t1,Q_RK(:,1)-M(:,2),'b')
title('Deviation link1')
legend('EULER','RK')

H5 = figure();
grid on
hold on
plot(t1,Q_EU(:,2)-M(:,3),'g')
plot(t1,Q_RK(:,2)-M(:,3),'b')
title('Deviation link2')
legend('EULER','RK')

% H6 = figure();
% grid on
% hold on
% plot(t1,Q_EU(:,3)-M(:,4),'g')
% plot(t1,Q_RK(:,3)-M(:,4),'b')
% title('Deviation link3')
% legend('EULER','RK')


end