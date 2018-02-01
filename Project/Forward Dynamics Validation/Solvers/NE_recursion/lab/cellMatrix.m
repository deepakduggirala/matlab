function [B,C] = cellMatrix()
clear
model = getModel(2);
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

a = zeros(3,NB);
b = zeros(3,NB);
c = zeros(3,NB);
d = zeros(3,NB);
e = zeros(1,NB);

C = zeros(26,1);
%calculating angular velocities
w(:,1)     = qd(1)*n(:,1);
for i=2:NB
    w(:,i) = qd(i)*n(:,i) + w(:,i-1);
end
%Calculating R.H.S of eqns
for i=1:NB
    a(:,i) = S(w(:,i))*(Icm{i}*w(:,i));
    b(:,i) = -m(i)*grav;
    c(:,i) = -S(w(:,i))*(qd(i)*n(:,i));
    if i~=NB
        d(:,i) = -S(w(:,i))*S(w(:,i))*r(:,i) - S(w(:,i+1))*S(w(:,i+1))*r(:,i+1);
    else
        d(:,i) = -S(w(:,i))*S(w(:,i))*r(:,i);
    end
    d(:,i) = d(:,i) - S(r(:,1))*c(:,1);
    e(i) = tau(i);
    a(:,i) = a(:,i) - S(r(:,i))*b(:,i);
    C((i-1)*13+1:i*13,1) = [a(:,i);b(:,i);c(:,i);d(:,i);e(i)];
end

A = zeros(13,13);
B = zeros(26,26);

for k=1:NB
    i=k;
    A(1:3,1:3) = eye(3);
    A(1:3,4:6) = -2*S(r(:,i));
    A(1:3,7:9) = Icm{i};
    A(1:3,10:12) = -m(i)*S(r(:,i));
    A(4:6,4:6) = eye(3);
    A(4:6,10:12) = m(i)*eye(3);
    A(7:9,7:9) = eye(3);
    A(7:9,13) = n(:,i);
    A(10:12,10:12) = eye(3);
    A(10:12,13) = -S(r(:,i))*n(:,i);
    A(13,1:3) = n(:,i)';
    B((k-1)*13+1:k*13,(k-1)*13+1:k*13) = A(:,:);
    if k<NB
        A = zeros(13,13);
        A(7:9,7:9) = -eye(3);
        A(10:12,7:9) = S(r(:,i)) + S(r(:,i+1));
        A(10:12,10:12) = -eye(3);
        B(1:13,14:26) = A(:,:);
        A = zeros(13,13);
        A(1:3,1:3) = -eye(3);
        A(4:6,4:6) = -eye(3);
        B(14:26,1:13) = A(:,:);
    end
    
end
end

