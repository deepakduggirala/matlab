function qdd_new = NE_rec4(model)
%Unpacking model
NB    = model.NB;
q     = model.q;
qd    = model.qd;
n0    = model.n0;
l0    = model.l0;
m     = model.mass;
Icm0  = model.Icm0;
grav     = -model.g;
tau   = model.tau;

%pre calculation of transformation matrices
R = preCalc(NB,n0,q);

%updating
l = zeros(3,NB);        
n = zeros(3,NB);
Xnet = eye(3);
Icm = cell(1,NB);
for i=1:NB
    n(:,i) = Xnet*n0(:,i);
    Xnet = Xnet*R(:,:,i);
    l(:,i) = Xnet*l0(:,i);
    Icm{i} = Xnet*Icm0(:,:,i)*Xnet';
end


%calculating angular velocities
w(:,1)     = qd(1)*n(:,1);
for i=2:NB
    w(:,i) = qd(i)*n(:,i) + w(:,i-1);
end

%flipping arrays
q = fliplr(q);
qd = fliplr(qd);
m = fliplr(m);
n = fliplr(n);
l = fliplr(l);
Icm = fliplr(Icm);
w = fliplr(w);
r = l/2;

a = zeros(3,NB);
b = zeros(3,NB);
c = zeros(3,NB);
d = zeros(3,NB);
e = zeros(1,NB);
%Calculating R.H.S of eqns
for i=1:NB
    a(:,i) = S(w(:,i))*(Icm{i}*w(:,i));
    b(:,i) = -m(i)*grav;
    c(:,i) = -S(w(:,i))*(qd(i)*n(:,i));
    if i~=NB
        d(:,i) = -S(w(:,i))*S(w(:,i))*r(:,i) - S(w(:,i+1))*S(w(:,i+1))*r(:,i+1);
    else
        d(:,i) = -S(w(:,i))*S(w(:,i))*r(:,i);
    end
    d(:,i) = d(:,i) - S(r(:,i))*c(:,i);
    e(i) = tau(i);
    a(:,i) = a(:,i) - S(r(:,i))*b(:,i);
end
K13 = cell(1,NB);
K14 = cell(1,NB);
K23 = cell(1,NB);
K24 = cell(1,NB);
D = cell(1,NB);
R53 = cell(1,NB);
R54 = cell(1,NB);
%base case n=1 for forward recursion
K13{1} = Icm{1};
K14{1} = -m(1)*S(r(:,1));
K23{1} = zeros(3);
K24{1} = m(1)*eye(3);
D{1} = n(:,1)'*(Icm{1} - m(1)*S(r(:,1))*S(r(:,1)))*n(:,1);
if NB~=1
    R53{1} = n(:,1)'*(-Icm{1} + m(1)*S(r(:,1))*(S(r(:,1)) + S(r(:,2))));
    R54{1} = -n(:,1)'*S(r(:,1))*m(1);
end

e(1) = e(1) - n(:,1)'*a(:,1) - 2*n(:,1)'*S(r(:,1))*b(:,1) + n(:,1)'*Icm{1}*c(:,1) + n(:,1)'*S(r(:,1))*m(1)*d(:,1);
%Forward recusrsion
for i=2:NB
    %Decoupling
    u_t = K13{i-1} + 2*S(r(:,i-1))*K23{i-1};
    v_t = K14{i-1} + 2*S(r(:,i-1))*K24{i-1};
    w_t = -u_t*n(:,i-1) + v_t*S(r(:,i-1))*n(:,i-1);
    K13{i} = Icm{i} + u_t - v_t*(S(r(:,i-1)) + S(r(:,i))) - w_t*R53{i-1}/D{i-1};
    K14{i} = -m(i)*S(r(:,i)) + v_t - w_t*R54{i-1}/D{i-1};
    a(:,i) = a(:,i) + a(:,i-1) + 2*S(r(:,i-1))*b(:,i-1) - u_t*c(:,i-1) -v_t*d(:,i-1) -w_t*e(i-1)/D{i-1};
    
    
    r_t = (-K23{i-1} + K24{i-1}*S(r(:,i-1)))*n(:,i-1)/D{i-1};
    K23{i} = K23{i-1} - K24{i-1}*(S(r(:,i-1)) + S(r(:,i))) - r_t*R53{i-1};
    K24{i} = m(i)*eye(3) + K24{i-1} - r_t*R54{i-1};
    
    b(:,i) = b(:,i) + b(:,i-1) -K23{i-1}*c(:,i-1) - K24{i-1}*d(:,i-1) - r_t*e(:,i-1);
    
    %Upper Triangularization
    q_t = -n(:,i)'*(K14{i} + 2*S(r(:,i))*K24{i});
    p_t = -n(:,i)'*(K13{i} + 2*S(r(:,i))*K23{i});
    D{i} = (-p_t + q_t*S(r(:,i)))*n(:,i);
    if i~=NB
        R53{i} = p_t - q_t*(S(r(:,i)) + S(r(:,i+1)));
        R54{i} = q_t;
    end
    e(i) = e(i) - n(:,i)'*a(:,i) - 2*n(:,i)'*S(r(:,i))*b(:,i) - p_t*c(:,i) - q_t*d(:,i);
    
end
% %backward recursion
qdd(NB) = e(NB)/D{NB};
acc(:,NB) = -(d(:,NB) + S(r(:,NB))*n(:,NB)*qdd(NB));
wd(:,NB) = -(c(:,NB) - n(:,NB)*qdd(NB));
for i=NB-1:-1:1
    qdd(i) = (e(i) + R53{i}*wd(:,i+1) + R54{i}*acc(:,i+1))/D{i};
    acc(:,i) = -(d(:,i) + S(r(:,i))*n(:,i)*qdd(i) + (S(r(:,i)) + S(r(:,i+1)))*wd(:,i+1) - acc(:,i+1));
    wd(:,i) = -(c(:,i) - wd(:,i+1) - n(:,i)*qdd(i));
end

qdd_new = fliplr(qdd);
qdd_new = qdd_new';
end


%skew symmetric matrix from vector to represent cross product
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
