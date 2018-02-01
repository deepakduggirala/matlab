function [tau, F, T, alpha, a, w, Icm, rc] = RNEA(model)

%Unpacking model
NB    = model.NB;
qdd   = model.qdd;
q     = model.q;
qd    = model.qd;
n0    = model.n0;
l0    = model.l0;
m     = model.mass;
F_ext = model.F_ext;
T_ext = model.T_ext;
Icm0  = model.Icm0;
g     = -model.g;


%initial conditions
w0    = [0 0 0]';
v0    = [0 0 0]';
alpha0= [0 0 0]';
a0    = [0 0 0]';


%pre calculation of transformation matrices
R = preCalc(NB,n0,q);

%updating
l = zeros(3,NB);        %homogeneous coords ?
n = zeros(3,NB);
Xnet = eye(3);
Icm = zeros(3,3,NB);
rc = zeros(3,NB);
for i=1:NB
    n(:,i) = Xnet*n0(:,i);
    Xnet = Xnet*R(:,:,i);
    l(:,i) = Xnet*l0(:,i);
    Icm(:,:,i) = Xnet*Icm0(:,:,i)*Xnet';
    rc(:,i) = l(:,i)/2;
end

w(:,1)     = w0 + qd(1)*n(:,1);
v(:,1)     = v0 + cross(w(:,1),l(:,1));
alpha(:,1) = alpha0 + qdd(1)*n(:,1) + qd(1)*cross(w0,n(:,1));
a(:,1)     = a0 + cross(alpha(:,1),l(:,1)) + cross(w(:,1),cross(w(:,1),l(:,1)));
a_cm(:,1)  = cross(alpha(:,1),rc(:,1)) + cross(w(:,1),cross(w(:,1),rc(:,1))) + a0;
for i=2:NB
    w(:,i) = qd(i)*n(:,i) + w(:,i-1);
    v(:,i) = v(:,i-1) + cross(w(:,i),l(:,i));
    
    alpha(:,i) = qdd(i)*n(:,i) + qd(i)*cross(w(:,i-1),n(:,i)) + alpha(:,i-1);
    a(:,i)     = a(:,i-1) + cross(alpha(:,i),l(:,i)) + cross(w(:,i),cross(w(:,i),l(:,i)));
    a_cm(:,i)  = cross(alpha(:,i),rc(:,i)) + cross(w(:,i),cross(w(:,i),rc(:,i))) + a(:,i-1);
end

% w,v,alpha,a,a_cm

F(:,NB) = m(NB)*a_cm(:,NB) + F_ext - m(NB)*g;
T(:,NB) = T_ext + Icm(:,:,NB)*alpha(:,NB) + cross(w(:,NB),Icm(:,:,NB)*w(:,NB)) - cross(-rc(:,NB),F(:,NB)) - cross(rc(:,NB),-F_ext);
for i=NB-1:-1:1
    F(:,i) = m(i)*a_cm(:,i) + F(:,i+1) - m(i)*g;
    T(:,i) = T(:,i+1) + Icm(:,:,i)*alpha(:,i) + cross(w(:,i),Icm(:,:,i)*w(:,i)) - cross(-rc(:,i),F(:,i)) - cross(rc(:,i),-F(:,i+1));
end


tau = zeros(NB,1);
for i=1:NB
    tau(i,1) = T(:,i)'*n(:,i);
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
    n(3),0,-n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end