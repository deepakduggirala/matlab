function [model1,model2] = halfParallelogram(tht)
%This function creates model which has the following fields of the
%parallelogram linkgae with revolute joints.

%model.length - model.length(i) is ith link length
%model.mass   - model.mass(i) is the ith link mass 
%model.q      - model.q(i) is joint angle parameter of ith link at start
%model.qd     - model.qd(i) is joint velocity parameter of ith link at start
%model.r0     - model.r0(:,i) is location of CoM in link coordinates
%model.n0     - model.n0(:,i) is unit vector along the axis of rotation at t=0;
%model.l0     - model.l0(:,i) is vector along ith link length
%model.theta  - model.theta is the angle link1 makes with the horiozontal
%at the start of the simulation

density = 7801;
model1.length = [0.1 0.1];
model2.length = [0.1 0.1];
model1.radius = [0.001 0.001];
model2.radius = [0.001 0.001];
tht1 = deg2rad(tht);
tht2 = tht1;
model1.theta = tht1;
model2.theta = tht2;

for i=1:2
    model1.mass(i) = density * model1.length(i)*pi*((model1.radius(i))^2);
end
for i=1:2
    model2.mass(i) = density * model2.length(i)*pi*((model2.radius(i))^2);
end

model1.qd= zeros(2,1);
model1.q = zeros(2,1);
model2.qd= zeros(2,1);
model2.q = zeros(2,1);

model1.n0 = [0 0
            0 0
            1 1];
model1.l0 = [model1.length(1)*cos(tht1)  model1.length(2)
             model1.length(1)*sin(tht1)  0
             0                          0];
        
model2.n0 = [0 0
            0 0
            1 1];
        
model2.l0 = [model2.length(1)*cos(tht2)  -model2.length(2)
             model2.length(1)*sin(tht2)  0
             0                          0];

model1.Icm0 = zeros(3,3,2);
model2.Icm0 = zeros(3,3,2);

model1.g     = [0 9.80665 0]';
model2.g     = [0 9.80665 0]';
model1.tau = zeros(2,1);
model2.tau = zeros(2,1);
model1.NB = 2;
model2.NB = 2;

R = zeros(3,3,2);
R(:,:,1) = rotMat(model1.n0(:,1), tht1);
R(:,:,2) = eye(3);

for i=1:2
    Ie = 0.5*model1.mass(i)*(model1.radius(i)^2);
    I  = model1.mass(i)*((model1.length(i)^2)/12) + Ie/2;
    J = [Ie 0 0
         0  I 0
         0  0 I];
    model1.Icm0(:,:,i) = R(:,:,i)*J*R(:,:,i)';
end

R = zeros(3,3,2);
R(:,:,1) = rotMat(model2.n0(:,1), tht2);
R(:,:,2) = eye(3);

for i=1:2
    Ie = 0.5*model2.mass(i)*(model2.radius(i)^2);
    I  = model2.mass(i)*((model2.length(i)^2)/12) + Ie/2;
    J = [Ie 0 0
         0  I 0
         0  0 I];
    model2.Icm0(:,:,i) = R(:,:,i)*J*R(:,:,i)';
end
end

function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end