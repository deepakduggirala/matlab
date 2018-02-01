dq1 = sym('dq1');
dq2 = sym('dq2');
dq3 = sym('dq3');

dx1 = sym('dx1');
dx2 = sym('dx2');
dx3 = sym('dx3');

m1 = sym('m1');
m2 = sym('m2');
m3 = sym('m3');

J = [dx1/dq1 dx1/dq2 dx1/dq3;dx2/dq1 dx2/dq2 dx2/dq3;dx3/dq1 dx3/dq2 dx3/dq3;];
Jt = [dx1/dq1 dx2/dq1 dx3/dq1;dx1/dq2 dx2/dq2 dx3/dq2;dx1/dq3 dx2/dq2 dx3/dq3;];

M=[m1 0 0;0 m2 0;0 0 m3];

x=M*J
Jt*x