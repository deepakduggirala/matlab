function [qdd,r] = mySolver(model)
l = model.l0(1,1);
g = -9.81;
qdd = 1.5*g*cos(model.q)/l;
r=zeros(2,1);
r(1,1)=l*cos(model.q)/2;
r(2,1) = l*sin(model.q)/2;
end