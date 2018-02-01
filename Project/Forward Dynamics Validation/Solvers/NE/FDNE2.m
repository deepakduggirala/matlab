function [qdd,A,B] = FDNE2(model)

%Unpacking model
NB    = model.NB;
q     = model.q;
qd    = model.qd;
n0    = model.n0;
l0    = model.l0;
m     = model.mass;
Icm0  = model.Icm0;
g     = -model.g;
tau     = model.tau;


%pre calculation of transformation matrices
R = preCalc(NB,n0,q);

%updating
l = zeros(3,NB);        
n = zeros(3,NB);
Xnet = eye(3);
Icm = zeros(3,3,NB);
IcmI = zeros(3,3,NB);
for i=1:NB
    n(:,i) = Xnet*n0(:,i);
    Xnet = Xnet*R(:,:,i);
    l(:,i) = Xnet*l0(:,i);
    Icm(:,:,i) = Xnet*Icm0(:,:,i)*Xnet';
    IcmI(:,:,i) = inv(Icm(:,:,i));
end
rc = l/2;


w(:,1)     = qd(1)*n(:,1);
for i=2:NB
    w(:,i) = qd(i)*n(:,i) + w(:,i-1);
end

A = zeros(6*NB,6*NB);
A(1:6,1:6) = getA(1,1,rc,l,IcmI,m);
for k = 2:NB
    ii = 6*(k-1)+1:6*k;
    jj = 6*(k-2)+1:6*(k-1);
    A(ii,ii) = getA(k,k,rc,l,IcmI,m);
    A(ii,jj) = getA(k,k-1,rc,l,IcmI,m);
    A(jj,ii) = getA(k-1,k,rc,l,IcmI,m);
end
B = zeros(6*NB,1);
for k=1:NB
    if(k~=NB)
        b = S(l(:,k))*IcmI(:,:,k)*(tau(k)*n(:,k) - tau(k+1)*n(:,k+1) - S(w(:,k))*Icm(:,:,k)*w(:,k));
    else
         b = S(l(:,k))*IcmI(:,:,k)*(tau(k)*n(:,k) - S(w(:,k))*Icm(:,:,k)*w(:,k));
    end
    B(6*(k-1)+1:6*(k-1)+3,1) = b-S(w(:,k))*S(w(:,k))*l(:,k);
    B(6*(k-1)+4:6*(k-1)+6,1) = m(k)*g;
end
% rcond(A)
x=A\B;
F =  zeros(3,NB);
for i=1:NB
    F(:,i) = x(6*(i-1)+4:6*i,1);
end

alpha = zeros(3,NB);
for i = 1:NB-1
    alpha(:,i) = IcmI(:,:,i)*(tau(i)*n(:,i) - tau(i+1)*n(:,i+1) - S(rc(:,i))*F(:,i) - S(rc(:,i))*F(:,i+1) - S(w(:,i))*(Icm(:,:,i)*w(:,i)));
end
alpha(:,NB) = IcmI(:,:,NB)*(tau(NB)*n(:,NB) - S(rc(:,NB))*F(:,NB) - S(w(:,NB))*(Icm(:,:,NB)*w(:,NB)));
% alpha
qdd = zeros(NB,1);
qdd(1,1) = n(:,1)'*(alpha(:,1))/((n(:,1)'*n(:,1)));
for i=2:NB
    qdd(i,1) = n(:,i)'*(alpha(:,i) - S(w(:,i-1))*(qd(i)*n(:,i)) - alpha(:,i-1))/((n(:,i)'*n(:,i)));
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
function X = S(n)
X = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
end
function X =getA(i,j,r,l,IcmI,m)
I = eye(3);
Z = zeros(3,3);
if(i == j)
    X = [-I S(l(:,i))*IcmI(:,:,i)*S(r(:,i));
        m(i)*I/2 -I];
elseif(i==j+1)
    X = [I Z;
        m(i)*I/2 Z];
elseif(i==j-1)
    X = [Z S(l(:,i))*IcmI(:,:,i)*S(r(:,i));
        Z I];
end
end