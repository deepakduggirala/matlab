function [thtdd,r] = EX2_solver(model)
    m1 = model.mass(1);
    m2 = model.mass(2);
    l = model.l0(1);
    tht1 = model.q(1);
    tht2 = model.q(2);
    thtd1 = model.qd(1);
    thtd2 = model.qd(2);
    tau = model.tau;
    g = 9.81;
    r0 = model.r0;
    NB = model.NB;
    
    R = zeros(4,4,2);
    R(:,:,1) = [cos(tht1) -sin(tht1) 0 0;
                sin(tht1) cos(tht1) 0  0;
                0 0 1 0;
                0 0 0 1];
    R(:,:,2) = [cos(tht2) -sin(tht2) 0 l;
                sin(tht2) cos(tht2) 0  0;
                0 0 1 0;
                0 0 0 1];
    
    r = zeros(4,NB);
    Xnet = eye(4);
    for i=1:NB
        Xnet = Xnet*R(:,:,i);
        r(:,i) = Xnet*r0(:,i);
    end
    
    D = [((m1*l^2)/3)+ (4*(m2*l^2)/3) + m2*cos(tht2)*l^2 ((m2*l^2)/3)+((m2*cos(tht2)*l^2)/2)
         ((m2*l^2)/3)+((m2*cos(tht2)*l^2)/2) (m2*l^2)/3];
    H = -0.5*m2*sin(tht2)*(l^2)*[thtd2^2 + 2*thtd1*thtd2; -thtd1^2];
    C = 0.5*g*l*[m1*cos(tht1)+m2*cos(tht1+tht2)+2*m2*cos(tht1); m2*cos(tht1+tht2)];
    
    thtdd = inv(D)*(tau - H -C);
    
end