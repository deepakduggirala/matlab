function [B,qdd] = openLoop(model)
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

r = l/2;

C = zeros(26,1);
%calculating angular velocities
w(:,1)     = qd(1)*n(:,1);
for i=2:NB
    w(:,i) = qd(i)*n(:,i) + w(:,i-1);
end


%Calculating R.H.S of eqns
for k=1:NB
    i=k;
    a = -S(w(:,i))*(Icm{i}*w(:,i));
    b = m(i)*grav;
    c = S(w(:,i))*(qd(i)*n(:,i));
    d = S(w(:,i))*S(w(:,i))*l(:,i);
    e = tau(i);
    C((k-1)*13+1:k*13,1) = [a;b;c;d;e];
end


B = zeros(13*NB,13*NB);

for k=1:NB
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
x = B\C;
qdd = zeros(NB,1);
for i=1:NB
    qdd(i) = x(13*i);
end
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
