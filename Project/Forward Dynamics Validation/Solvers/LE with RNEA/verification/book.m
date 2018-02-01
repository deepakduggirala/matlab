function book(q,qd,qdd)
n1 = [0 1 0]';
n2 = [sin(q(2)) 0 cos(q(2))]';

w2 = qd(2)*n1
w32 = qd(3)*n2;
w3 = w2+w32

a(:,1)=zeros(3,1);
alpha(:,1) = zeros(3,1);


A2 = [cos(q(2)) 0 sin(q(2))
        0         1     0       
      -sin(q(2)) 0   cos(q(2))];

A32 = [cos(q(3)) -sin(q(3)) 0;
       sin(q(3)) cos(q(3)) 0;
        0         0         1];
A3 = A2*A32;

UQ2b = [0 -0.05 0]';
UQ2 = A2*UQ2b;

UP2b = [0 0.05 0]';
UP2 = A2*UP2b;

UP3b = [0 -0.05 0]';
UP3 = A3*UP3b;

alpha(:,2) = [0 1 0]'*qdd(2);
a(:,2) = S(UQ2)*n1*qdd(2) -S(w2)*(S(w2)*UQ2);


g3 = S(w2)*(S(w2)*UP2) -S(w3)*(S(w3)*UP3);
gtht3 = qd(2)*qd(3)*[cos(q(2)) 0 -sin(q(2))]';
alpha(:,3) = alpha(:,2) + qdd(3)*n2 + gtht3;
a(:,3) = a(:,2) + (S(UP3) - S(UP2))*alpha(:,2) + S(UP3)*n2*qdd(3) + g3+S(UP3)*gtht3;

a,alpha
end

function X = S(n)
X = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
end