function [Qdd,r,w,RHS] = FDNEQ(model)
NB    = model.NB;
Q     = model.Q;
Qd    = model.Qd;
l0    = model.l0;
m     = model.mass;
Icm0  = model.Icm0;
grav  = model.g;

% n0    = model.n0;
% q     = model.q;

H = [0 0 0
     1 0 0
     0 1 0
     0 0 1];

%updating with quaternions
Q_net = [1;0;0;0];
l = zeros(3,NB);
Icm = cell(1,NB);
for i=1:NB
    Q_net = quat_prod(Q_net,Q(:,i));
    l(:,i) = Q_on_v(Q_net,l0(:,i));
    R = quat2RotMat(Q_net);
    Icm{i} = R*Icm0(:,:,i)*R';
end

% %updating with rotation matrices
% 
% %pre calculation of transformation matrices
% R = preCalc(NB,n0,q);
% 
% l_R = zeros(3,NB);        
% n = zeros(3,NB);
% Xnet = eye(3);
% Icm_R = cell(1,NB);
% for i=1:NB
%     n(:,i) = Xnet*n0(:,i);
%     Xnet = Xnet*R(:,:,i);
%     l_R(:,i) = Xnet*l0(:,i);
%     Icm_R{i} = Xnet*Icm0(:,:,i)*Xnet';
% end
% 
% %testing
% r1 = 1;
% r2 = 1;
% tol = 10^(-15);
% x = abs(l - l_R);
% err_l = max(reshape(x,[],1));
% if(err_l > tol)
%     r1 = 0;
% end
% for i=1:NB
%     x = abs(Icm{i} - Icm_R{i});
%     err_I = max(reshape(x,[],1));
%     if(err_I > tol)
%         r2 = 0;
%         break;
%     end
% end

r = l/2;

RHS = zeros(13*NB,1);
LHS = zeros(13*NB,13*NB);

%calculating angular velocities
w = zeros(3,NB);
w(:,1) = H'*2*quat_prod(Qd(:,1),qConj(Q(:,1)));
for i = 2:NB
    w(:,i) = w(:,i-1) + H'*2*quat_prod(Qd(:,i),qConj(Q(:,i)));
end

%Calculating R.H.S of eqns
for i=1:NB
    a = -S(w(:,i))*(Icm{i}*w(:,i));
    b = m(i)*grav;
    c = S(w(:,i))*S(w(:,i))*l(:,i);
    if i == 1
        d = -T(quat_prod(Qd(:,i),qConj(Q(:,i))))*H*(w(:,i));
    else
        d = -T(quat_prod(Qd(:,i),qConj(Q(:,i))))*H*(w(:,i) - w(:,i-1));
    end
    RHS((i-1)*13+1:i*13,1) = [a;b;c;d]; 
end

for i=1:NB
    A = zeros(13,13);
    A(1:3,1:3) = Icm{i};
    A(1:3,4:6) = S(r(:,i));
    A(4:6,4:6) = -eye(3);
    A(4:6,7:9) = (m(i)/2)*eye(3);
    A(7:9,1:3) = S(l(:,i));
    A(7:9,7:9) = eye(3);
    A(10:13,1:3) = H;
    A(10:13,10:13) = -2*T(qConj(Q(:,i)));
    LHS((i-1)*13+1:i*13,(i-1)*13+1:i*13) = A(:,:);
    if i<NB
        A = zeros(13,13);
        A(1:3,4:6) = S(r(:,i));
        A(4:6,4:6) = eye(3);
        LHS((i-1)*13+1:i*13,(i)*13+1:(i+1)*13) = A(:,:);
        
        A = zeros(13,13);
        A(4:6,7:9) = (m(i)/2)*eye(3);
        A(7:9,7:9) = -eye(3);
        A(10:13,1:3) = -H;
        LHS((i)*13+1:(i+1)*13,(i-1)*13+1:i*13) = A(:,:);
    end
end
Qdd = zeros(4,NB);
x = LHS\RHS;
for i=1:NB
    Qdd(:,i) = x(13*(i-1)+10:13*(i-1)+13);
end
end

function x = quat_prod(a,b)     %Quaternion Conjugate
s1 = a(1);
s2 = b(1);
v1 = a(2:4);
v2 = b(2:4);
x = [s1*s2 - dot(v1,v2); s1*v2 + s2*v1 + cross(v1,v2)];
end

function x = Q_on_v(Q,v)        %Applying quaternion Q on vector v: QvQ-1
x = quat_prod(Q,[0;v]);
x = quat_prod(x,qConj(Q));
x = x(2:4);
end

function Q1 = qConj(Q)          %Quaternion conjugate
Q1 = [Q(1);-Q(2);-Q(3);-Q(4)];
end

function R = quat2RotMat(q)     %Rotation matrix corresponding to a quaternion Q
    R = [q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2, 2*(q(2)*q(3) - q(1)*q(4)),         2*(q(2)*q(4) + q(1)*q(3))
         2*(q(2)*q(3) + q(1)*q(4))          q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2, 2*(q(3)*q(4) - q(1)*q(2))
         2*(q(2)*q(4) - q(1)*q(3)),         2*(q(3)*q(4) + q(1)*q(2)),         q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2];
end

function X = T(q)
    s = q(1);
    v = q(2:4);
    X = [s, -v'; v, s*eye(3) - S(v)];
end

function X = S(n)
    X = [0,-n(3),n(2);
         n(3),0,-n(1);
        -n(2),n(1),0];
end
    
function R = preCalc(n,n0,q)
R = zeros(3,3,n);
for i = 1:n
    R(:,:,i) = rotMat(n0(:,i),q(i));
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