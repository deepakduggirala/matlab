function model = helix()
r=1;
c=1;

t = 0:pi/30:3*pi;
x = r*sin(t);
y = r*cos(t);
z = c*t;
N = size(t);
N = N(2);
model.NB = N-1;
l0 = zeros(3,N-1);
for i=1:N-1
    src = [x(i),y(i),z(i)]';
    dst = [x(i+1),y(i+1),z(i+1)]';
    l0(:,i) = dst-src;
end
model.l0 = l0;

n0 = zeros(3,N-1);
n0(2,1) = 1;
for i=2:N-1
    n = cross(l0(:,i),l0(:,i-1));
    n = n/norm(n);
    n0(:,i) = n;
end
model.n0 = n0;

density = 7801;
for i = 1:model.NB
    model.length(i) = 0.1*norm(l0(:,i));
end
model.radius = model.length/100;
for i=1:model.NB
    model.mass(i) = density * model.length(i)*pi*((model.radius(i))^2);
end

model.qd= zeros(model.NB,1);
model.q = zeros(model.NB,1);

J = zeros(3,3,model.NB);
Ie = 0.5*model.mass(1)*(model.radius(1)^2);
I  = model.mass(1)*((model.length(1)^2)/12) + Ie/2;
J1        = [Ie 0 0
             0 I 0
             0 0 I];
l_x = [1,0,0]';
for i=1:model.NB
    tht = acos(dot(l_x, l0(:,i))/norm(l0(:,i)));
    n = cross(l_x,l0(:,i));
    n = n/norm(n);
    R = rotMat(n,tht);
    J(:,:,i) = R*J1*R';
end
model.Icm0 = J;

% model.g     = [0 0 9.80665]';
model.g = zeros(3,1);
model.tau = zeros(model.NB,1);
model.F_ext = zeros(model.NB,1);
model.T_ext = zeros(model.NB,1);
model.K_t = 0.1;
end

function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end