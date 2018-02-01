function qdd = CCFDNE_planar(model)
NB    = model.NB;       %no:of links (scalar) : 1x1
q     = model.q;        %Joint angle (scalar) : NBx1
qd    = model.qd;       %Joint velocity (scalar) : NBx1
n0    = model.n0;       %Joint axis vector (3x1) : 3xNB
l0    = model.l0;       %length vectors(3x1) : 3xNB
m     = model.mass;     %Mass (scalar) : NBx1
Icm0  = model.Icm0;     %Moment of Inertia Tensor(3x3) : 3x3xNB
grav  = -model.g;       %acceleration due to graviry (vector) : 3x1
tau   = model.tau;      %Joint torque (scalar) : NBx1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRE-PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[l0,Icm, grav] = preProcess(l0,Icm0,n0,NB, grav);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pre calculation of transformation matrices
R = preCalc(NB,q);

%updating
l = zeros(2,NB);        
Xnet = eye(2);
for i=1:NB
    Xnet = Xnet*R(:,:,i);
    l(:,i) = Xnet*l0(1:2,i);
end

r = l/2;

w = zeros(NB,1);
w(1)     = qd(1);
for i=2:NB
    w(i) = qd(i) + w(i-1);
end

B = zeros(7*NB,1);
for k=1:NB
    i=k;
    a = 0;
    b = m(i)*grav(1:2);
    c = 0;
    d = [-l(1,i)*w(i)*w(i);-l(2,i)*w(i)*w(i)];
    e = tau(i);
    B((k-1)*7+1:k*7,1)=[a;b;c;d;e];
end
B = [B;zeros(3,1);tau(NB+1)];


A = zeros(7*NB, 7*NB);
for k=1:NB
    i=k;
    T = zeros(7,7);
    T(1,1) = -1;
    T(1,2) = -r(2,i);
    T(1,3) = r(1,i);
    T(1,4) = Icm(3,3,i);
    T(2,2) = -1;
    T(2,5) = m(i)/2;
    T(3,3) = -1;
    T(3,6) = m(i)/2;
    T(4,4) = 1;
    T(4,7) = -1;
    T(5,4) = l(2,i);
    T(5,5) = 1;
    T(6,4) = -l(1,i);
    T(6,6) = 1;
    T(7,1) = 1;
    A((k-1)*7+1:k*7,(k-1)*7+1:k*7) = T(:,:);
     if k<NB
        T = zeros(7,7);
        T(1,1) = 1;
        T(1,2) = -r(2,i);
        T(1,3) = r(1,i);
        T(2,2) = 1;
        T(3,3) = 1;
        A((k-1)*7+1:k*7,(k)*7+1:(k+1)*7) = T(:,:);
        T = zeros(7,7);
        T(2,5) = m(i+1)/2;
        T(3,6) = m(i+1)/2;
        T(4,4) = -1;
        T(5,5) = -1;
        T(6,6) = -1;                                                                                                                                                                                                              1;
        A((k)*7+1:(k+1)*7,(k-1)*7+1:k*7) = T(:,:);
     end
end

C = zeros(4,7*NB);
C(1,(NB-1)*7+4) = 1;
C(2,(NB-1)*7+5) = 1;
C(3,(NB-1)*7+6) = 1;

D = zeros(7*NB,4);
D(7*(NB-1)+1,1) = 1;
D(7*(NB-1)+1,2) = -r(2,NB);
D(7*(NB-1)+1,3) = r(1,NB);
D(7*(NB-1)+2,2) = 1;
D(7*(NB-1)+3,3) = 1;

E = zeros(4,4);
E(1,4) = 1;
E(4,1) = 1;

A_c = [A D;C E];
X = A_c\B;

qdd = zeros(NB+1,1);
for i=1:NB
    qdd(i) = X(7*i);
end
qdd(NB+1) = X(7*NB+4);
end


function R = preCalc(n,q)
R = zeros(2,2,n);
for i = 1:n
    R(:,:,i) = [cos(q(i)) -sin(q(i)); sin(q(i)) cos(q(i))];
end
end

function [l,Icm,grav] = preProcess(l0,Icm0,n0,NB,grav)
k = n0(:,1)/norm(n0(:,1));
e = l0(:,1)/norm(l0(:,1));
t = cross(k,e);

RX = [e t k]';

l = zeros(3,NB);
Icm = zeros(3,3,NB);
for i=1:NB
    l(:,i) = RX*l0(:,i);
    Icm(:,:,i) = RX*Icm0(:,:,i)*RX';
end
grav = RX*grav;
end