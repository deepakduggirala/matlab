function model = getModel(NB)
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

%constructing I
Ie = 0.5*model.mass(1)*(model.radius(1)^2);
I  = model.mass(i)*((model.length(i)^2)/12) + Ie/2;
Icm      = [Ie 0 0
             0 I 0
             0 0 I];
c = [0.05 0 0];
I = cell(1,NB);
for i=1:NB
    I{i} = mcI(model.mass(i),c,Icm);
end
model.I = I;


%Constructing coordinate transformations XT
%Link to link ctransformations at t=0
Xtree = cell(1,NB);
Xtree{1} = eye(6);
% t = [0 0 1
%      0 -1 0
%      1 0 0];
t = eye(3);
temp = [t' zeros(3)
       zeros(3) t'];
Xtree{2} = temp*xlt([0.1 0 0]);
% t = [0 -1 0
%     0 0 -1
%     1 0 0];
t = eye(3);
temp = [t' zeros(3)
       zeros(3) t'];
Xtree{3} = temp*xlt([0.1 0 0]);
model.Xtree = Xtree;

%gravity
model.g     = [0;0;0;0;-9.80665;0];

%joint constraint
model.s = [0 0 1 0 0 0]';
model.qd= zeros(NB,1);
model.q = zeros(NB,1);
model.tau = {0,0,0};
end