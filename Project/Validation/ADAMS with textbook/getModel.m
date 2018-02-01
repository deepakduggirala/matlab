function model = getModel( NB )
%This function creates model which has the following fields of the
%linkgae(unbranched) with revolute joints.

%model.NB     - No:Of links
%model.length - model.length(i) is ith link length
%model.mass   - model.mass(i) is the ith link mass 
%model.tau    - model.tau(i) is the torque on ith link about axis n = R*n0;
%model.q      - model.q(i) is joint angle parameter of ith link at start
%model.qd     - model.qd(i) is joint velocity parameter of ith link at start
%model.r0     - model.r0(:,i) is location of CoM in link coordinates
%model.n0     - model.n0(:,i) is unit vector along the axis of rotation at t=0;
%model.l0     - model.l0(:,i) is vector along ith link length
%model.J      - model.J(:,:,i) is the Inertia matrix of ith link


model.NB = NB;

linearDensity = 7801;
model.length = 0.1*ones(NB,1);
model.mass = linearDensity * model.length;

model.tau = zeros(NB,1);
model.qd= zeros(NB,1);
model.q = zeros(NB,1);

r0 = zeros(4,NB);
r0(1,:) = 0.5*model.length;
r0(4,:) = 1;
model.r0 = r0;

n0 = zeros(4,NB);
n0(3,:) = 1;
model.n0 = n0;

l0 = zeros(3,NB);
l0(1,:) = model.length;           %TODO for any angle
model.l0 = l0;

J0 = zeros(4,4,NB);
for i=1:NB
    j = zeros(4,4);
    Iyy = model.mass(i)*((model.length(i)^2)/3);
    Ixx = 0;
    Izz = model.mass(i)*((model.length(i)^2)/3);
    j(1,1) = (Iyy+Izz-Ixx)/2;
    j(2,2) = (Izz+Ixx-Iyy)/2;
    j(3,3) = (Ixx+Iyy-Izz)/2;
    j(:,4) = 0;
    j(4,:) = 0;
    j(4,4) = model.mass(i);
    j(1,4) = model.mass(i)*model.length(i)/2;
    j(4,1) = model.mass(i)*model.length(i)/2;
    J0(:,:,i) = j;          %J is a 3d matrix
end
model.J=J0;
end

