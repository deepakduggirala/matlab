function model = getModelN(NB)
%This function creates model which has the following fields of the
%linkgae(unbranched) with revolute joints.

%model.NB     - No:of links
%model.length - model.length(i) is ith link length
%model.mass   - model.mass(i) is the ith link mass 
%model.q      - model.q(i) is joint angle parameter of ith link at start
%model.qd     - model.qd(i) is joint velocity parameter of ith link at start
%model.r0     - model.r0(:,i) is location of CoM in link coordinates
%model.n0     - model.n0(:,i) is unit vector along the axis of rotation at t=0;
%model.l0     - model.l0(:,i) is vector along ith link length
%model.T_ext  - model.T_ext is an external torque on the manipulator(end of the last link)
%model.F_ext  - model.F_ext is an external force on the manipulator


model.NB = NB;
density = 7801;
model.length = 0.1*ones(NB,1);
model.radius = model.length/100;
for i=1:NB
    model.mass(i) = density * model.length(i)*pi*((model.radius(i))^2);
end

model.qd= zeros(NB,1);
model.q = zeros(NB,1);


l0 = zeros(3,NB);
l0(1,:) = 0.1;
model.l0=l0;

n0 = zeros(3,NB);
n0(3,:) = 1;
model.n0=n0;
model.r0=l0/2;

J = zeros(4,4,NB);
j = zeros(4,4);
i=1;
Ie = model.mass(i)*(model.radius(i)^2)/2;
I  = model.mass(i)*((model.length(i)^2)/3) + Ie/2;
    j(1,1) = I;
    j(4,4) = model.mass(i);
    j(1,4) = model.mass(i)*model.length(i)/2;
    j(4,1) = model.mass(i)*model.length(i)/2;
for i = 1:NB
    J(:,:,i) = j;
end
model.Icm0 = J;
model.tau = zeros(NB,1);
end