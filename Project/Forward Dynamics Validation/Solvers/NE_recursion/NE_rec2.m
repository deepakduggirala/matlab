function qdd = NE_rec2(model)
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
% Icm = zeros(3,3,NB);
for i=1:NB
    n(:,i) = Xnet*n0(:,i);
    Xnet = Xnet*R(:,:,i);
    l(:,i) = Xnet*l0(:,i);
    Icm{i} = Xnet*Icm0(:,:,i)*Xnet';
end
r = l/2;
%calculating angular velocities
w(:,1)     = qd(1)*n(:,1);
for i=2:NB
    w(:,i) = qd(i)*n(:,i) + w(:,i-1);
end

%R.H.S matrix before gaussian elimination
a = zeros(3,NB);
b = zeros(3,NB);
c = zeros(3,NB);
d = zeros(3,NB);
e = zeros(1,NB);
g = zeros(NB,3);
h = zeros(NB,3);
H = zeros(NB);
RR21 = zeros(3,3,NB);
for i=1:NB
    a(:,i) = m(i)*grav;
    b(:,i) = S(w(:,i))*S(w(:,i))*l(:,i);
    c(:,i) = S(w(:,i))*(qd(i)*n(:,i));
    d(:,i) = -S(w(:,i))*(Icm{i}*w(:,i));
    e(i) = tau(i);
end

%Forward Elimination
d(:,1) = d(:,1)+S(r(:,1))*a(:,1) - m(i)*S(r(:,1))*b(:,1)/2;
e(1) = e(1) + n(:,1)'*d(:,1);
H(1) = (n(:,1)')*(Icm{1}-m(1)*S(r(:,1))*S(r(:,1)))*n(:,1);
g(1,:) = 2*n(:,1)'*S(r(:,1));
h(1,:) = n(:,1)';
K44(:,:,1) = -eye(3);
K33(:,:,1) = eye(3);
K22(:,:,1) = eye(3);
K11(:,:,1) = -eye(3);
K14(:,:,1) = zeros(3);
K24(:,:,1) = zeros(3);
K34(:,:,1) = zeros(3);
K45(:,:,1) = (Icm{1}-m(1)*S(r(:,1))*S(r(:,1)))*n(:,1);
RR31(:,:,1) = zeros(3);
RR41(:,:,1) = 2*S(r(:,1));
for i=2:NB
    L13 = -m(i)*inv(K22(:,:,i-1))*S(l(:,i-1))/2;
    L14 = -m(i)*inv(K22(:,:,i-1))*K24(:,:,i-1)/2;
    L23 = inv(K22(:,:,i-1))*S(l(:,i-1));
    L24 = inv(K22(:,:,i-1))*K24(:,:,i-1);
    L34 = inv(K33(:,:,i-1))*K34(:,:,i-1);
    L35 = -inv(K33(:,:,i-1))*n(:,i-1);
    L14 = L14 - L13*inv(K33(:,:,i-1))*K34(:,:,i-1);
    L15 = L13*inv(K33(:,:,i-1))*n(:,i-1);
    L24 = L24 - L23*inv(K33(:,:,i-1))*K34(:,:,i-1);
    L25 = L23*inv(K33(:,:,i-1))*n(:,i-1);
    L15 = L15 - L14*inv(K44(:,:,i-1))*K45(:,:,i-1);
    L25 = L25 - L24*inv(K44(:,:,i-1))*K45(:,:,i-1);
    k31 = inv(K33(:,:,i-1))*RR31(:,:,i-1);
    k21 = -L23*inv(K33(:,:,i-1))*RR31(:,:,i-1);
    k11 = -eye(3) - L13*inv(K33(:,:,i-1))*RR31(:,:,i-1);
    k11 = k11 - L14*inv(K44(:,:,i-1))*RR41(:,:,i-1);
    k11 = k11 - L15*inv(H(i-1))*g(i-1,:);
    k14 = -L14*inv(K44(:,:,i-1)) - L15*inv(H(i-1))*h(i-1,:);
    k21 = k21 - L24*inv(K44(:,:,i-1))*RR41(:,:,i-1);
    k21 = k21 - L25*inv(H(i-1))*g(i-1,:);
    k24 = -L24*inv(K44(:,:,i-1)) - L25*inv(H(i-1))*h(i-1,:);
    L35 = L35 - L34*inv(K44(:,:,i-1))*K45(:,:,i-1);
    k31 = k31 - L34*inv(K44(:,:,i-1))*RR41(:,:,i-1);
    k31 = k31 - L35*inv(H(i-1))*g(i-1,:);
    k34 = -L34*inv(K44(:,:,i-1)) - L35*inv(H(i-1))*h(i-1,:);
    
    
    k22 = eye(3) - k21*inv(k11)*m(i)/2;
    k24 =  k24 - k21*inv(k11)*k14;
    R21 = -k21*inv(k11);
    k32 = -k31*inv(k11)*m(i)/2;
    k33 = eye(3) - k32*inv(k22)*S(l(:,i));
    k34 = k34 -k31*inv(k11)*k14 - k32*inv(k22)*k24;
    R31 = -k31*inv(k11) - k32*inv(k22)*R21;
    k42 = -S(r(:,i))*inv(k11)*m(i)/2;
    k44 = -eye(3) - S(r(:,i))*inv(k11)*k14;
    R41 = S(r(:,i)) - S(r(:,i))*inv(k11);
    k43 = Icm{i} - k42*inv(k22)*S(l(:,i));
    k44 = k44 - k42*inv(k22)*k24;
    R41 = R41 - k42*inv(k22)*R21;
    k44 = k44 - k43*inv(k33)*k34;
    k45 = k43*inv(k33)*n(:,i);
    R41 = R41 - k43*inv(k33)*R31;
    H(i)   = -n(:,i)'*inv(k44)*k45;
    g(i,:) = -n(:,i)'*inv(k44)*R41;
    h(i,:) = -n(:,i)'*inv(k44);
    
    a(:,i) = a(:,i) - m(i)*inv(K22(:,:,i-1))*b(:,i-1)/2 -L13*inv(K33(:,:,i-1))*c(:,i-1);
    a(:,i) = a(:,i) - L14*inv(K44(:,:,i-1))*d(:,i-1) - L15*inv(H(i-1))*e(i-1);
    b(:,i) = b(:,i) + inv(K22(:,:,i-1))*b(:,i-1) -L23*inv(K33(:,:,i-1))*c(:,i-1);
    b(:,i) = b(:,i) - L24*inv(K44(:,:,i-1))*d(:,i-1) - L25*inv(H(i-1))*e(i-1);
    c(:,i) = c(:,i) + inv(K33(:,:,i-1))*c(:,i-1) - L34*inv(K44(:,:,i-1))*d(:,i-1) - L35*inv(H(i-1))*e(i-1);
    
    b(:,i) = b(:,i) - k21*inv(k11)*a(:,i);
    c(:,i) = c(:,i) - k31*inv(k11)*a(:,i) - k32*inv(k22)*b(:,i);
    d(:,i) = d(:,i) - S(r(:,i))*inv(k11)*a(:,i) - k42*inv(k22)*b(:,i) - k43*inv(k33)*c(:,i);
    e(:,i) = e(:,i) - n(:,i)'*inv(k44)*d(:,i);
    
    K44(:,:,i) = k44;
    K33(:,:,i) = k33;
    K22(:,:,i) = k22;
    K11(:,:,i) = k11;
    K14(:,:,i) = k14;
    K24(:,:,i) = k24;
    K34(:,:,i) = k34;
    K45(:,:,i) = k45;
    RR21(:,:,i) = R21;
    RR31(:,:,i) = R31;
    RR41(:,:,i) = R41;
end


%Backward substituition
qdd(NB)   = inv(H(NB))*(e(NB));
T(:,NB)   = inv(K44(:,:,NB))*(d(:,NB) - K45(:,:,NB)*qdd(NB));
wd(:,NB)  = inv(K33(:,:,NB))*(c(:,NB) - K34(:,:,NB)*T(:,NB) + n(:,NB)*qdd(NB));
acc(:,NB) = inv(K22(:,:,NB))*(b(:,NB) - S(l(:,NB))*wd(:,NB)-K24(:,:,NB)*T(:,NB));
F(:,NB)   = inv(K11(:,:,NB))*(a(:,NB) - m(NB)*acc(:,NB)/2 - K14(:,:,NB)*T(:,NB));
for i = NB-1:-1:1
    qdd(i)   = inv(H(i))*(e(i) - g(i,:)*F(:,i+1) - h(i,:)*T(:,i+1));
    T(:,i)   = inv(K44(:,:,i))*(d(:,i) - K45(:,:,i)*qdd(i) - RR41(:,:,i)*F(:,i+1) - T(:,i+1));
    wd(:,i)  = inv(K33(:,:,i))*(c(:,i) - K34(:,:,i)*T(:,i) + n(:,i)*qdd(i) - RR31(:,:,i)*F(:,i+1));
    acc(:,i) = inv(K22(:,:,i))*(b(:,i) - S(l(:,i))*wd(:,i)-K24(:,:,i)*T(:,i) - RR21(:,:,i)*F(:,i+1));
    F(:,i)   = inv(K11(:,:,i))*(a(:,i) - m(i)*acc(:,i)/2 - K14(:,:,i)*T(:,i) - F(:,i+1)); %#ok<*MINV>
end
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

function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end
