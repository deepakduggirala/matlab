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
%model.k_t    - model.K_t is the spring constant of torsional springs at the joints
%model.K_d    - model.K_d is the damping coefficients of the torsional damper at the joints 


model.NB = NB;
density = 7801;
model.length = 0.1*ones(NB,1);
model.radius = model.length/100;
for i=1:NB
    model.mass(i) = density * model.length(i)*pi*((model.radius(i))^2);
end

model.qd= zeros(NB,1);
model.q = zeros(NB,1);

%3link
if(NB==3)
    l0 = [0.1 0.1 0.1;
           0  0   0;
           0  0   0];
    n0 = [0 0 0;
          0 0 0;
          1 1 1];
end

%2link
if (NB ==2)
    l0 = [0.1 0.1;
           0  0;
           0  0];
    n0 = [0 0;
          0 0;
          1 1];    
end

%single link
if(NB==1)
    l0 = [0.1 0 0]';
    n0 = [0 0 1]';
end

model.l0=l0;
model.n0=n0;
model.r0 = l0/2;
J = zeros(3,3,NB);
Ie = 0.5*model.mass(1)*(model.radius(1)^2);
I  = model.mass(1)*((model.length(1)^2)/12) + Ie/2;
J(:,:,1) = [Ie 0 0
             0 I 0
             0 0 I];
J(:,:,2) = [Ie 0 0
            0 I 0
            0 0 I];
J(:,:,3) = [Ie 0 0
             0 I 0
             0 0 I];

model.Icm0 = J;

model.g     = [0 9.80665 0]';
model.tau = zeros(NB,1);
model.F_ext = zeros(3,1);
model.T_ext = zeros(3,1);
model.K_t = 100;
model.K_d = 10;
end