function [q,qd] = motionIntegration(q,qd,qdd,timeStep)
    qd = qd + qdd*timeStep;
    q = q + qd*timeStep;
end