function B = CCFDNE(model)
NB    = model.NB;
q     = model.q;
qd    = model.qd;
n0    = model.n0;
l0    = model.l0;
m     = model.mass;
Icm0  = model.Icm0;
grav     = -model.g;
tau   = model.tau;

%pre calculation of transformation matrices
R = preCalc(NB,n0,q);

%updating
l = zeros(3,NB);        
n = zeros(3,NB);
Xnet = eye(3);
Icm = cell(1,NB);
for i=1:NB
    n(:,i) = Xnet*n0(:,i);
    Xnet = Xnet*R(:,:,i);
    l(:,i) = Xnet*l0(:,i);
    Icm{i} = Xnet*Icm0(:,:,i)*Xnet';
end

%calculating angular velocities
w(:,1)     = qd(1)*n(:,1);
for i=2:NB
    w(:,i) = qd(i)*n(:,i) + w(:,i-1);
end

r = l/2;

C = zeros(13*NB+7,1);

%Calculating R.H.S of eqns
for k=1:NB-1
    i=k;
    a = -S(w(:,i))*(Icm{i}*w(:,i));
    b = m(i)*grav;
    c = S(w(:,i))*(qd(i)*n(:,i));
    d = S(w(:,i))*S(w(:,i))*l(:,i);
    e = tau(i);
    C((k-1)*13+1:k*13,1) = [a;b;c;d;e];
end

a = -2*S(w(:,NB))*Icm{NB}*w(:,NB) - 2*S(r(:,NB))*m(NB)*grav;
b = 2*m(NB)*grav;
c = S(w(:,NB-1))*qd(NB)*n(:,NB);
d = S(w(:,NB))*S(w(:,NB))*l(:,NB);
e = tau(NB);
f = zeros(3,1);
g = -S(w(:,NB))*S(w(:,NB))*l(:,NB);
h = tau(NB+1);

C(13*(NB-1)+1:13*NB+7,1) = [a;b;c;d;e;f;g;h];

B = zeros(13*NB+7,13*NB+7);

for k=1:NB-1
    i=k;
    A = zeros(13,13);
    A(1:3,1:3) = -eye(3);
    A(1:3,4:6) = S(r(:,i));
    A(1:3,7:9) = Icm{i};
    A(4:6,4:6) = -eye(3);
    A(4:6,10:12) = m(i)*eye(3)/2;
    A(7:9,7:9) = eye(3);
    A(7:9,13) = -n(:,i);
    A(10:12,7:9) = S(l(:,i));
    A(10:12,10:12) = eye(3);
    A(13,1:3) = n(:,i)';
    B((k-1)*13+1:k*13,(k-1)*13+1:k*13) = A(:,:);
    if k<NB
        A = zeros(13,13);
        A(1:3,1:3) = eye(3);
        A(1:3,4:6) = S(r(:,i));
        A(4:6,4:6) = eye(3);
        B((k-1)*13+1:k*13,(k)*13+1:(k+1)*13) = A(:,:);
        A = zeros(13,13);
        A(4:6,10:12) = m(i+1)*eye(3)/2;
        A(7:9,7:9) = -eye(3);
        A(10:12,10:12) = -eye(3);
        B((k)*13+1:(k+1)*13,(k-1)*13+1:k*13) = A(:,:);
    end
end
k = NB-1;
B(k*13+1:k*13+3,(k-1)*13+10:(k-1)*13+12) = S(r(:,NB))*m(NB);

A = zeros(13,13);
A(1:3,1:3) = -eye(3);
A(1:3,4:6) = S(r(:,NB));
A(1:3,7:9) = 2*Icm{NB};
A(1:3,10:12) = S(r(:,NB))*m(NB);
A(4:6,4:6) = -eye(3);
A(4:6,10:12) = m(NB)*eye(3);
A(7:9,7:9) = eye(3);
A(7:9,13) = -n(:,NB);
A(10:12,7:9) = S(l(:,NB));
A(10:12,10:12) = eye(3);
A(13,1:3) = n(:,NB)';
B((NB-1)*13+1:NB*13,(NB-1)*13+1:NB*13) = A(:,:);

A = zeros(13,7);
A(1:3,1:3) = -eye(3);
A(1:3,4:6) = -S(r(:,NB));
A(4:6,4:6) = -eye(3);
B((NB-1)*13+1:NB*13,(NB)*13+1:(NB)*13+7) = A(:,:);

A = zeros(7,13);
A(1:3,7:9) = eye(3);
% A(4:6,7:9) = -S(l(:,NB));
A(4:6,10:12) = eye(3);
B((NB)*13+1:(NB)*13+7,(NB-1)*13+1:NB*13) = A(:,:);

A = zeros(7,7);
A(1:3,7) = -n0(:,NB+1);
A(7,1:3) = n0(:,NB+1)';
B((NB)*13+1:(NB)*13+7,(NB)*13+1:(NB)*13+7) = A(:,:);

B(13*NB+4:13*NB+6,(NB-2)*13+10:(NB-2)*13+12) = -eye(3)/2;

B;
% x = B\C
% qdd = zeros(NB,1);
% for i=1:NB
%     qdd(i) = x(13*i);
% end
end

%skew symmetric matrix from vector to represent cross product
function X = S(n)
X = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
end

function R = preCalc(n,n0,q)
R = zeros(3,3,n);
for i = 1:n
    R(:,:,i) = rotMat(n0(:,i),q(i));
end
end

%Rodriguez rotation matrix
function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end
