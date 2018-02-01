clear;
clc;
c1 = sym('cos(q1)');
c2 = sym('cos(q2)');
s1 = sym('sin(q1)');
s2 = sym('sin(q2)');

R1 = [c1 -s1 0;s1 c1 0;0 0 1];
R2 = [c2 -s2 0;s2 c2 0;0 0 1];

l1l = sym('l');
% l2l = sym('l2');
l10 = [l1l;0;0];
l20 = [l1l;0;0];
l1 = R1*l10;
l2 = R1*R2*l20;

qd1 = sym('qd1');
qd2 = sym('qd2');
qdd1 = sym('qdd1');
qdd2 = sym('qdd2');

n1 = [0;0;1];
n2 = [0;0;1];

w1 = qd1*n1;
w2 = w1 + qd2*n2;

v1 = cross(w1,l1);
v2 = v1 + cross(w2,l2);

wd1 = qdd1*n1;
wd2 = qdd2*n2 + qd2*cross(w1,n2) + wd1;

a1 = -cross(l1,wd1) + cross(w1,cross(w1,l1));
a2 = a1 - cross(l2,wd2) + cross(w2,cross(w2,l2));

a1cm = a1/2;
a2cm = (a1+a2)/2;

% m2 = sym('m2');
m1 = sym('m');
m2=m1;
g = sym('g');

F2 = m2*a2cm - m2*g;
F1 = m1*a1cm + F2 -m1*g;

I1 = m1*(l1l*l1l)/12;
I2 = m2*(l1l*l1l)/12;

Icm1 = [0 0 0;0 I1 0;0 0 I1];
Icm2 = [0 0 0;0 I2 0;0 0 I2];

T2 = Icm2*wd2 + cross(w2,Icm2*w2) + cross(l2/2,F2);
T1 = Icm1*wd1 + cross(w1,Icm1*w1) + cross(l1/2,F1) + cross(l1/2,F2) + T2;