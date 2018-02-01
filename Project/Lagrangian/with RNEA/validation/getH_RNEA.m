function tau = getH_RNEA(model_LE)
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

NB = model_LE.NB;
model.NB = model_LE.NB;

model.length = model_LE.length;  %m
model.mass = model_LE.mass;%kg

model.qdd = zeros(NB,1);
model.qd= model_LE.qd;
model.q = model_LE.q;

model.n0 = model_LE.n0;
model.l0 = model_LE.l0;

J = zeros(3,3,NB);
for i=1:NB
    j = zeros(3,3);
    Iyy = model.mass(i)*((model.length(i)^2)/3);
    Ixx = 0;
    Izz = model.mass(i)*((model.length(i)^2)/3);
    
    j(1,1) = Ixx;
    j(2,2) = Iyy;
    j(3,3) = Izz;
J(:,:,i) = j;    
end
model.Icm0 = J;

model.F_ext = zeros(3,1);
model.T_ext = zeros(3,1);
model.g     = [0 0 0]';


tau = RNEA_bd(model);
end