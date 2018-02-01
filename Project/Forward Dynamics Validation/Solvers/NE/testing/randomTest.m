clc
clear all;
NB=2;
q = (pi/2)*rand(NB,1);
qd = 10*rand(NB,1);
% q=zeros(NB,1);
% qd=zeros(NB,1);
qdd = 10*rand(NB,1);

model = getModel(NB);
model.q = q;
model.qd = qd;
model.qdd = qdd;

[tau, F, T, alpha, a, w, Icm, rc] = RNEA(model);
[A, B,w2,Icm2] = FDNE(model);

x = zeros(9*NB,1);
for k=1:NB
    x(9*(k-1)+1:9*(k-1)+3,1) = a(:,k);
    x(9*(k-1)+4:9*(k-1)+6,1) = F(:,k);
    x(9*(k-1)+7:9*(k-1)+9,1) = alpha(:,k);
end
y=A*x;
x1 = -cross(w(:,1),Icm(:,:,1)*w(:,1));
x2 = -cross(w(:,2),Icm(:,:,2)*w(:,2));
y1 = cross(rc(:,1),F(:,1)) + cross(rc(:,1),F(:,2)) + Icm(:,:,1)*alpha(:,1);
y2 = cross(rc(:,2),F(:,2)) + Icm(:,:,2)*alpha(:,2);


y(7:9)  %A*x(7:9)
y1      %RNEA A*x
B(7:9)  %B(7:9)
x1      %RNEA B
y(16:18)
y2
B(16:18)
x2

abs(y-B)<10^-6