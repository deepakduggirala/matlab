function [qdd,r] = LE_solver(model)
%solves the forward dynamics problem with Lagrangian-Euler method

%unpacking model
NB = model.NB;
q = model.q;
qd = model.qd;
J = model.J;
r0 = model.r0;
l0 = model.l0;
n0 = model.n0;
m = model.mass;
tau = model.tau;

%Pre-calculation and storing of these terms avoids recalculating.
%U(:,:,i,j) is equivalent to getU(i,j)
%R(:,:,i) is equivalent to getR(i)
%RD(:,:,i) is equivalent to getRD(i)
[R,RD,U] = preCalc(NB,l0,n0,q);

r = zeros(4,NB);
Xnet = eye(4);
for i=1:NB
    Xnet = Xnet*R(:,:,i);
    r(:,i) = Xnet*r0(:,i);
end

D = getD(NB,U,J);
H = getH(NB,qd,J,U,R,RD,n0,q);
C = getC(NB,r0,m,U);
qdd = inv(D)*(tau-H+C);
end

function [R,RD,U] = preCalc(n,l0,n0,q)
R = zeros(4,4,n);
RD = zeros(4,4,n);
for i = 1:n
    R(:,:,i) = getR(i,l0,n0,q);
    RD(:,:,i) = getRD(i,n0,q);
end

U = zeros(4,4,n,n);
for i=1:n
    for j = 1:i
        U(:,:,i,j) = getU(i,j,R,RD);
    end
end
end
function T = getU(i,j,R,RD)
T = eye(4);
for ii = 1:i
    if ii==j
        T = T*RD(:,:,ii);
    else
        T=T*R(:,:,ii);
    end
end
end
function R = rotMatDD(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,n(1);
    -n(2),n(1),0];
R = (n*n')*cos(tht) - eye(3)*cos(tht) - S*sin(tht); 
end
function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end
function R = rotMatD(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,n(1);
    -n(2),n(1),0];
R = (n)*(n')*(sin(tht)) - eye(3)*(sin(tht)) + (S)*cos(tht);
end
function T = getR(i,l0,n0,tht)
T = zeros(4,4);
T(1:3,1:3) = rotMat(n0(1:3,i),tht(i));
T(4,4) = 1;
if i == 1
    T(1:3,4) = [0 0 0]';
else
    T(1:3,4) = l0(:,i-1);
end
end
function T = getRD(i,n0,tht)
T = zeros(4);
T(1:3,1:3) = rotMatD(n0(1:3,i),tht(i));
end
function T = getRDD(i,n0,tht)
T = zeros(4);
T(1:3,1:3) = rotMatDD(n0(1:3,i),tht(i));
end
function R = getU2(i,j,k,X,XD,n0,tht)
R = eye(4);
if i > max(j,k)
    if j==k
        for ii=1:i
            if ii==j
                R = R*getRDD(ii,n0,tht);
            else
                R = R*X(:,:,ii);
            end
        end
    else
        for ii=1:i
            if ii==j
                R = R*XD(:,:,ii);
            end
            if ii==k
                R=R*XD(:,:,ii);
            end
            if ii ~= j && ii ~= k
                R = R*X(:,:,ii);
            end
        end
    end
else
    if i == j && i ~= k
        for ii=1:i
            if ii==j
                R = R*getRDD(ii,n0,tht);
            else
                R = R*X(:,:,ii);
            end
        end
    end
    if i == k && i~= j
        for ii=1:i
            if ii==k
                R = R*getRDD(ii,n0,tht);
            else
                R = R*X(:,:,ii);
            end
        end
    end
    if i == j && i== k
        for ii=1:i
            if ii==j
                R = R*getRDD(ii,n0,tht);
            else
                R = R*X(:,:,ii);
            end
        end
    end
end
end

function D = getD(n,U,J)
%Calculates the acceleration related Inertia matrix("D").
%this function exploits the fact that D is a symmetrical matrix by
%calculating only the lower triangle terms and for upper triangle terms are
%copied from thier counterparts in lower triangle.
%
D = zeros(n);
    for i = 1:n
        for k = 1:i
            sum = 0;
            for j = i:n
                sum = sum + trace(U(:,:,j,k)*J(:,:,j)*(U(:,:,j,i)'));
            end
            D(i,k) = sum;
        end
    end
    
    
%Upper triangular terms
    for i=1:n-1
        for j = i+1:n
            D(i,j) = D(j,i);
        end
    end
end

function H = getH(nL,thtd,J,U,R,RD,n0,tht)
H = zeros(nL,1);
for i=1:nL
    sum = 0;
    for k = 1:nL
        for m = 1:nL
            sum = sum + get_h(nL,i,k,m,J,U,R,RD,n0,tht)*thtd(k)*thtd(m);
        end
    end
    H(i) = sum;
end
end
function sum =get_h(nL,i,k,m,J,U,R,RD,n0,tht)
sum = 0;
for j=max(i,max(k,m)):nL
    sum = sum + trace(getU2(j,k,m,R,RD,n0,tht)*J(:,:,j)*(U(:,:,j,i)'));
end
end

function C = getC(nL,r0,m,U)
C = zeros(nL,1);
for i=1:nL
    sum = 0;
    for j=i:nL
      sum = sum + m(j)*[0 -9.8 0 0]*U(:,:,j,i)*r0(:,j);
    end
    C(i) = sum;
end
end