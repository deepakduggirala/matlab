% clear
clc
NB=2;
q = (pi/2)*rand(NB,1);
qd = 10*rand(NB,1)

m1 = getModel1(NB);
m  = getModel(NB);

m.q=q;
m.qd=qd;

m1.q=q;
m1.qd=qd;

[qdd1,D1,H1,C1] = LE_solver(m1);
[qdd,D,H,C] = FDLE(m);
% qdd-qdd1
% D-D1
H,H1
% C+C1