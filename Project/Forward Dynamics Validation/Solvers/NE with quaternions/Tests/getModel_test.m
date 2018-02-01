function model = getModel_test(NB)
%This function creates model which has the following fields of the
%linkgae(unbranched) with spherical joints.

%model.NB     - No:of links
%model.length - model.length is link length
%model.mass   - model.mass is link mass 
%model.Q      - model.Q(:,i) is rotation quaternion of ith link at start
%model.Qd     - model.qd(i) is derivative of rotation quaternion of ith link at start
%model.r0     - model.r0(:,i) is location of CoM in link coordinates
%model.l0     - model.l0(:,i) is vector along ith link length
%model.T_ext  - model.T_ext is an external torque on the manipulator(end of the last link)
%model.F_ext  - model.F_ext is an external force on the manipulator


model.NB = NB;
density = 7801;
model.length = 0.1;
model.radius = model.length/100;
model.g      = [0 -9.80665 0]';
model.F_ext  = zeros(3,1);
model.T_ext  = zeros(3,1);
model.mass = ones(NB,1)*(density * model.length*pi*((model.radius)^2));

model.Qd = zeros(4,NB);
model.Q = zeros(4,NB);

model.n0 = [0 1 0
            0 0 0
            1 0 1];
model.q  = [deg2rad(30); deg2rad(30); deg2rad(30)];

for i = 1:NB
    model.Q(:,i) = quat_from_angle_axis(model.q(i),model.n0(:,i));
end

%3link
if(NB==3)
    l0 = [0.1 0     0;
           0  0     0.1;
           0  0.1   0];
end

%2link
if (NB ==2)
    l0 = [0.1 0.1;
           0  0;
           0  0];   
end

%single link
if(NB==1)
    l0 = [0.1 0 0]';
end

model.l0=l0;
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
end

function Q = quat_from_angle_axis(tht,n)
    Q = [cos(tht/2);sin(tht/2)*[n(1); n(2); n(3)]];
end