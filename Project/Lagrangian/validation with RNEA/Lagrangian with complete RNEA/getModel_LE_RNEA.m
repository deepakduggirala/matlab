function model = getModel_LE_RNEA(NB)
model.NB = NB;

linearDensity = 1;
model.length = 0.1*ones(NB,1);
model.mass = linearDensity * model.length;

model.qd= zeros(NB,1);
model.q = zeros(NB,1);
model.qdd = zeros(NB,1);

model.tau = zeros(NB,1);

% if(NB==3)
%     l0=[0.1 0 0.1;
%         0   0 0;
%         0 0.1 0];
%     n0=[0 1 0;
%         0 0 1;
%         1 0 0];
%     r0=[0.05 0 0.05;
%         0    0    0;
%         0 0.05    0];
% end
% 
% if(NB == 2)
%     l0=[0.1 0;
%         0   0;
%         0 0.1];
%     n0=[0 0;
%         0 1;
%         1 0];
%     r0=[0.05 0;
%         0    0;
%         0 0.05];
% end
% 
% r0(4,:) = 1;
% model.l0=l0;
% model.r0=r0;
% model.n0=n0;

r0 = zeros(4,NB);
r0(1,:) = 0.5*model.length;
r0(4,:) = 1;
model.r0 = r0;

n0 = zeros(3,NB);
n0(3,:) = 1;
model.n0 = n0;

l0 = zeros(3,NB);
l0(1,:) = model.length;           %TODO for any angle
model.l0 = l0;

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