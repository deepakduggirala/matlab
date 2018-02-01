function [thtdd,r] = EX_solver(model)

m = model.mass;
r = model.r0;
tht = model.q;
thtd = model.qd;
r1 = 0.05;
r2 = 0.05;
l1 = 0.1;
l2 = 0.1;
IZ1 = m(1)*(l1^2)/3;
IZ2 = m(2)*(l2^2)/3;
tau = model.tau;

alpha = IZ1 + IZ2 + m(1)*r1^2 + m(2)*(l1^2 + r2^2);
beta = m(2)*l1*r2;
del = IZ2 + m(2)*r2^2;

D = [alpha + 2*beta*cos(tht(2)) del+beta*cos(tht(2));del+beta*cos(tht(2)) del];
H = [-beta*sin(tht(2))*thtd(2) -beta*sin(tht(2))*(thtd(1)+thtd(2)); beta*sin(tht(2))*thtd(1) 0];
thtdd = inv(D)*(tau - H*thtd);
    
end