function FDNE_rev_sph(model)
NB = model.NB;
l0    = model.l0;
Icm0  = model.Icm0;
n0    = model.n0;

%update l, Icm
l = zeros(3,NB);
n = zeros(3,NB);
Icm = cell(1,NB);
Xnet = eye(3);
for i=1:NB
    n(:,i) = Xnet*n0(:,i);
    Xnet = Xnet*getTransformation(model, i);
    l(:,i) = Xnet*l0(:,i);
    Icm{i} = Xnet*Icm0(:,:,i)*Xnet';
end

r = l/2;

RHS = zeros(13*NB,1);
LHS = zeros(13*NB,13*NB);

%update angular velocities
w = zeros(3,NB);
w(:,1) = getRelAngVel(model, n, 1);
for i=2:NB
    w(:,i) = w(:,i-1) + getRelAngVel(model, n, i);
end

%Calculating R.H.S of eqns
for i=1:NB
    RHS((i-1)*13+1:i*13,1) = getRHS(model, Icm, l, w, n, i);
end

%Calculating L.H.S of eqns
for i=1:NB
    LHS((i-1)*13+1:i*13,(i-1)*13+1:i*13) = getDiagBlock(i);
    
    
end
end

function X = getTransformation(model, i)
if(strcmpi(model.joint{i},'rev'))
    X = rotMat(model.n0(:,i),model.tht(i));
elseif(strmcpi(model.joint{i},'sph'))
    X = quat2RotMat(model.Q(:,i));
end
end

%Rodriguez rotation matrix
function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end

function R = quat2RotMat(q)     %Rotation matrix corresponding to a quaternion Q
    R = [q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2, 2*(q(2)*q(3) - q(1)*q(4)),         2*(q(2)*q(4) + q(1)*q(3))
         2*(q(2)*q(3) + q(1)*q(4))          q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2, 2*(q(3)*q(4) - q(1)*q(2))
         2*(q(2)*q(4) - q(1)*q(3)),         2*(q(3)*q(4) + q(1)*q(2)),         q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2];
end

function w_rel = getRelAngVel(model, n, i)
H = [0 0 0
     1 0 0
     0 1 0
     0 0 1];
 
if(strcmpi(model.joint{i},'rev'))
    w_rel = model.thtd(i)*n(:,i);
elseif(strmcpi(model.joint{i},'sph'))
    w_rel = H'*2*quat_prod(model.Qd(:,1),qConj(model.Q(:,1)));
end
end

function x = quat_prod(a,b)     %Quaternion Conjugate
s1 = a(1);
s2 = b(1);
v1 = a(2:4);
v2 = b(2:4);
x = [s1*s2 - dot(v1,v2); s1*v2 + s2*v1 + cross(v1,v2)];
end

function Q1 = qConj(Q)          %Quaternion conjugate
Q1 = [Q(1);-Q(2);-Q(3);-Q(4)];
end

function X = S(n)
    X = [0,-n(3),n(2);
         n(3),0,-n(1);
        -n(2),n(1),0];
end

function RHS = getRHS(model, Icm, l, w, n, i)
if(strcmpi(model.joint{i},'rev'))
    a = -S(w(:,i))*(Icm{i}*w(:,i));
    b = model.mass(i)*model.grav;
    c = S(w(:,i))*(model.thtd(i)*n(:,i));
    d = S(w(:,i))*S(w(:,i))*l(:,i);
    e = model.tau(i);
    RHS = [a;b;c;d;e];
elseif(strmcpi(model.joint{i},'sph'))
    a = -S(w(:,i))*(Icm{i}*w(:,i));
    b = model.mass(i)*model.grav;
    c = S(w(:,i))*S(w(:,i))*l(:,i);
    if i == 1
        d = -T(quat_prod(model.Qd(:,i),qConj(model.Q(:,i))))*H*(w(:,i));
    else
        d = -T(quat_prod(model.Qd(:,i),qConj(model.Q(:,i))))*H*(w(:,i) - w(:,i-1));
    end
    RHS((i-1)*13+1:i*13,1) = [a;b;c;d];
end
end