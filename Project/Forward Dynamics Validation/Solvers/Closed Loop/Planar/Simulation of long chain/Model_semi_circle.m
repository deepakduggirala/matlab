function model = Model_semi_circle(NB, R)
%NB links whose endpoints are arranged along a semi circle of radius R in XY
%plane
%R     - radius of the ring (semi circle) on which the links are arranged 
%model.NB     - Number of Bodies or Links
%model.mass   - model.mass(i) is the ith link mass 
%model.q      - model.q(i) is joint angle parameter of ith link at start
%model.qd     - model.qd(i) is joint velocity parameter of ith link at start
%model.n0     - model.n0(:,i) is unit vector along the axis of rotation at t=0;
%model.l0     - model.l0(:,i) is vector along ith link length
%model.Icm0   - model.Icm0(:,:,i) is a 3x3 matrix of moment of Inertia
%model.g   - Acceleration due to gravity vctor(column)

model.NB = NB;
model.q = zeros(NB);
model.qd = zeros(NB);
model.g     = [0 9.80665 0]';
model.tau = zeros(NB+1);

density = 7801;
del_tht = pi/NB;
model.l0 = zeros(2,NB);
model.n0 = zeros(3,NB+1);
model.Icm0 = zeros(3,3,NB);

A = [R;0];
tht = 0;
for i = 1:NB
    B = R*[cos(tht+del_tht);sin(tht+del_tht)];
    model.l0(:,i) = B - A;
    A = B;
    tht = tht + del_tht;
end

%plot_ring_linkage(model,R)

for i=1:NB+1
    model.n0(:,i) = [0;0;1];
end

length = norm(model.l0(:,1));
radius = length/100;
link_mass = density*length*pi*(radius^2);
Ie = 0.5*link_mass*(radius^2);
I = link_mass*((length^2)/12) + Ie/2;
J = [Ie 0 0
     0  I 0
     0  0 I];
X_axis = [1;0];
for i=1:NB
    model.mass(i) = link_mass;
    tht = acos(dot(X_axis,model.l0(:,i)/length));
    R_mat = rotMat(model.n0(:,i),tht);
    model.Icm0(:,:,i) = R_mat*J*R_mat';
end
end

function plot_ring_linkage(model,R)
hold on
A = [R;0];
for i = 1:model.NB
    B = A + model.l0(:,i);
    plot([A(1,1),B(1,1)],[A(2,1),B(2,1)]);
    A = B;
end
end

function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end