function model = Parallelogram(tht)
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
model.length = [0.1 0.2 0.1];
model.radius = [0.001 0.001 0.001];
tht = deg2rad(tht);
model.theta = tht;

for i=1:3
    model.mass(i) = density * model.length(i)*pi*((model.radius(i))^2);
end

model.qd= zeros(3,1);
model.q = zeros(3,1);

model.n0 = [0 0 0 0
            0 0 0 0
            1 1 1 1];
model.l0 = [model.length(1)*cos(tht) model.length(2) -model.length(3)*cos(tht)
            model.length(1)*sin(tht) 0 -model.length(3)*sin(tht)
            0 0 0];

R = zeros(3,3,3);
R(:,:,1) = rotMat(model.n0(:,1), tht);
R(:,:,2) = eye(3);
R(:,:,3) = rotMat(model.n0(:,3),mod(pi+tht,2*pi));

model.Icm0 = zeros(3,3,3);

model.g     = [0 9.80665 0]';
model.tau = zeros(4,1);
model.NB = 3;


for i=1:3
    Ie = 0.5*model.mass(i)*(model.radius(i)^2);
    I  = model.mass(i)*((model.length(i)^2)/12) + Ie/2;
    J = [Ie 0 0
         0  I 0
         0  0 I];
    model.Icm0(:,:,i) = R(:,:,i)*J*R(:,:,i)';
end
end

function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end