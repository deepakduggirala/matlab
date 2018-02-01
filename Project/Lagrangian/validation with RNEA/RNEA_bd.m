function [tau] = RNEA_bd(model)

%Unpacking model
NB    = model.NB;
qdd   = model.qdd;
q     = model.q;
qd    = model.qd;
n0    = model.n0;
l0    = model.l0;
m     = model.mass;
Icm0  = model.Icm0;
g     = model.g;


z0 = [0 0 1]';
rc = l0/2;

%initial conditions
w0    = [0 0 0]';
alpha0= [0 0 0]';
a0    = g;


%pre calculation of transformation matrices
R = preCalc(NB,n0,q);



w(:,1) = R(:,:,1)'*(qd(1)*z0 + w0);
alpha(:,1) = R(:,:,1)'*(qdd(1)*z0 + qd(1)*cross(w0,z0) + alpha0);
a(:,1)     = R(:,:,1)'*a0 + cross(alpha(:,1),l0(:,1)) + cross(w(:,1),cross(w(:,1),l0(:,1)));
a_cm(:,1)  = cross(alpha(:,1),-rc(:,1)) + cross(w(:,1),cross(w(:,1),-rc(:,1))) + a(:,1);
for i=2:NB
    w(:,i) = R(:,:,i)'*(qd(i)*z0 + w(:,i-1));
    alpha(:,i) = R(:,:,i)'*(qdd(i)*z0 + qd(i)*cross(w(:,i-1),z0) + alpha(:,i-1));
    a(:,i)     = R(:,:,i)'*a(:,i-1) + cross(alpha(:,i),l0(:,i)) + cross(w(:,i),cross(w(:,i),l0(:,i)));
    a_cm(:,i)  = cross(alpha(:,i),-rc(:,i)) + cross(w(:,i),cross(w(:,i),-rc(:,i))) + a(:,i);
end

i=NB;
f(:,i) = m(i)*a_cm(:,i);
n(:,i) = cross(rc(:,i),m(i)*a_cm(:,i)) + Icm0(:,:,i)*alpha(:,i) + cross(w(:,i),Icm0(:,:,i)*w(:,i));
for i=NB-1:-1:1
    f(:,i) = R(:,:,i+1)*(f(:,i+1)) + m(i)*a_cm(:,i);
    n(:,i) = R(:,:,i+1)*(n(:,i+1) + cross(R(:,:,i+1)'*l0(:,i),f(:,i+1))) + cross(rc(:,i),m(i)*a_cm(:,i)) + Icm0(:,:,i)*alpha(:,i) + cross(w(:,i),Icm0(:,:,i)*w(:,i));
end

tau = zeros(NB,1);
for i=1:NB
    tau(i,1) = n(:,i)'*R(:,:,i)'*z0;
end

end

function R = preCalc(n,n0,q)
R = zeros(3,3,n);
for i = 1:n
    R(:,:,i) = rotMat(n0(:,i),q(i));
end
end
function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end
