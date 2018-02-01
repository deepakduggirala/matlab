function [A,B,w,Icm] = FDNE(model)

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
for i=1:NB
    n(:,i) = Xnet*n0(:,i);
    Xnet = Xnet*R(:,:,i);
    l(:,i) = Xnet*l0(:,i);
    Icm(:,:,i) = Xnet*Icm0(:,:,i)*Xnet';
end
rc = l/2;


w(:,1)     = qd(1)*n(:,1);
for i=2:NB
    w(:,i) = qd(i)*n(:,i) + w(:,i-1);
end

A = zeros(9*NB,9*NB);
A(1:9,1:9) = getA(1,1,rc,l,Icm,m);
for k = 2:NB
    ii = 9*(k-1)+1:9*k;
    jj = 9*(k-2)+1:9*(k-1);
    A(ii,ii) = getA(k,k,rc,l,Icm,m);
    A(ii,jj) = getA(k,k-1,rc,l,Icm,m);
    A(jj,ii) = getA(k-1,k,rc,l,Icm,m);
end
B = zeros(9*NB,1);
for k=1:NB
    B(9*(k-1)+1:9*(k-1)+3,1) = S(w(:,k))*S(w(:,k))*l(:,k);
    B(9*(k-1)+4:9*(k-1)+6,1) = m(k)*g;
    if(k~=NB)
        B(9*(k-1)+7:9*(k-1)+9,1) = tau(k)*n(:,k) - tau(k+1)*n(:,k+1) - S(w(:,k))*Icm(:,:,k)*w(:,k);
    else
        B(9*(k-1)+7:9*(k-1)+9,1) = tau(k)*n(:,k) - S(w(:,k))*Icm(:,:,k)*w(:,k);
    end
end
% rcond(A)
% x=A\B;
% alpha = zeros(3,NB);
% for i = 1:NB
%     alpha(:,i) = x(9*(i-1)+7:9*i,1);
% end
% % alpha
% qdd = zeros(NB,1);
% qdd(1,1) = n(:,1)'*(alpha(:,1))/((n(:,1)'*n(:,1)));
% for i=2:NB
%     qdd(i,1) = n(:,i)'*(alpha(:,i) - S(w(:,i-1))*(qd(i)*n(:,i)) - alpha(:,i-1))/((n(:,i)'*n(:,i)));
% end
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
function X =getA(i,j,r,l,Icm,m)
I = eye(3);
Z = zeros(3,3);
if(i == j)
    X = [I Z S(l(:,i));
        m(i)*I/2 -I Z;
        Z S(r(:,i)) Icm(:,:,i)];
elseif(i==j+1)
    X = [-I Z Z;
        m(i)*I/2 Z Z;
        Z Z Z];
elseif(i==j-1)
    X = [Z Z Z;
        Z I Z;
        Z S(r(:,i)) Z];
end
end