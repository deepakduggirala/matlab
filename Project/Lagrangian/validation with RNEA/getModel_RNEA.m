function model = getModel_RNEA(NB)
%This function creates model which has the following fields of the
%linkgae(unbranched) with revolute joints.

%model.NB     - No:Of links
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

linearDensity = 1;
model.length = 0.1*ones(NB,1);
model.mass = linearDensity * model.length;

model.qdd = zeros(NB,1);
model.qd= zeros(NB,1);
model.q = zeros(NB,1);

% r0 = zeros(3,NB);
% r0(1,:) = 0.5*model.length;
% model.r0 = r0;
% 
% n0 = zeros(3,NB);
% n0(3,:) = 1;
% model.n0 = n0;
% 
% l0 = zeros(3,NB);
% l0(1,:) = model.length;           %TODO for any angle
% model.l0 = l0;

% l0 = 0.1*[1 0;0 0;0 1];
% n0 = [0 0;0 1;1 0];

% model.l0=l0;
% model.n0=n0;



J = zeros(3,3,NB);
for i=1:NB
    j = zeros(3,3);
    Iyy = model.mass(i)*((model.length(i)^2)/12);
    Ixx = 0;
    Izz = model.mass(i)*((model.length(i)^2)/12);
    
    j(1,1) = Ixx;
    j(2,2) = Iyy;
    j(3,3) = Izz;
J(:,:,i) = j;    
end
model.Icm0 = J;

model.F_ext = zeros(3,1);
model.T_ext = zeros(3,1);
model.g     = [0 9.80665 0]';
end



