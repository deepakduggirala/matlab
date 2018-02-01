for t=0:delT:timeStep
    m    = m_RK;
    q1   = m.q;
    qd1  = m.qd;
    qdd1 = LE_solver(m);
    
    q2   = q1  + 0.5*qd1*delT;
    qd2  = qd1 + 0.5*qdd1*delT;
    m.q  = q2;
    m.qd = qd2;
    qdd2 = LE_solver(m);
    
    q3   = q1  + 0.5*qd2*delT;
    qd3  = qd1 + 0.5*qdd2*delT;
    m.q  = q3;
    m.qd = qd3;
    qdd3 = LE_solver(m);
    
    q4   = q1  + qd3*delT;
    qd4  = qd1 + qdd3*delT;
    m.q  = q4;
    m.qd = qd4;
    qdd4 = LE_solver(m);
    
    qf   = q1  + (delT/6)*(qd1 + 2*qd2 + 2*qd3 + qd4);
    qdf  = qd1 + (delT/6)*(qdd1 + 2*qdd2 + 2*qdd3 + qdd4);
    m_RK.q = qf;
    m_RK.qd = qdf;
end