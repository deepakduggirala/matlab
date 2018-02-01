%author: Deepak Duggirala
%Date:   27-08-2015
%Description: O(n) Recursive Forward Dynamics with Newton -Euler Equations

function qdd = solver()
model = getModelSin(100);
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
r = l/2;
%calculating angular velocities
w(:,1)     = qd(1)*n(:,1);
for i=2:NB
    w(:,i) = qd(i)*n(:,i) + w(:,i-1);
end

a = zeros(3,NB);
b = zeros(3,NB);
c = zeros(3,NB);
d = zeros(3,NB);
e = zeros(1,NB);
K33{i} = cell(1,NB);
K34{i} = cell(1,NB);
K35{i} = cell(1,NB);
R31{i} = cell(1,NB);
R32{i} = cell(1,NB);
K44{i} = cell(1,NB);
K45{i} = cell(1,NB);
R41{i} = cell(1,NB);
R42{i} = cell(1,NB);
K55{i} = cell(1,NB);
R51{i} = cell(1,NB);
R52{i} = cell(1,NB);
%Calculating R.H.S of eqns
for i=1:NB
    a(:,i) = -S(w(:,i))*(Icm{i}*w(:,i));
    b(:,i) = m(i)*grav;
    c(:,i) = S(w(:,i))*(qd(i)*n(:,i));
    d(:,i) = S(w(:,i))*S(w(:,i))*l(:,i);
    e(i) = tau(i);
end
%base case n=1 for forward recursion
K33{1} = eye(3);
K34{1} = zeros(3);
K35{1} = -n(:,1);
R31{1} = zeros(3);
R32{1} = zeros(3);
K44{1} = eye(3);
K45{1} = S(l(:,1))*n(:,1);
R41{1} = zeros(3);
R42{1} = zeros(3);
K55{1} = n(:,1)'*Icm{1}*n(:,1) - n(:,1)'*S(r(:,1))*m(i)*S(r(:,1))*n(:,1);
R51{1} = n(:,1)';
R52{1} = 2*n(:,1)'*S(r(:,1));
d(:,1) = d(:,1) - S(l(:,1))*c(:,1);
e(:,1) = e(:,1) + n(:,1)'*a(:,1) + n(:,1)'*S(r(:,1))*b(:,1) - n(:,1)'*Icm{1}*c(:,1) - n(:,1)'*S(r(:,1))*m(1)*d(:,1)/2;

for i=2:NB
    %decoupling blocks
    L34 = K33{i-1}\K34{i-1};
    L35 = K33{i-1}\K35{i-1};
    k31 = K33{i-1}\R31{i-1};
    k32 = K33{i-1}\R32{i-1};
    L35 = L35 - L34*(K44{i-1}\K45{i-1});
    k31 = k31 - L34*(K44{i-1}\R41{i-1});
    k32 = k32 - L34*(K44{i-1}\R42{i-1});
    L45 = K44{i-1}\K45{i-1};
    k41 = K44{i-1}\R41{i-1};
    k42 = K44{i-1}\R42{i-1};
    k31 = k31 - L35*(K55{i-1}\R51{i-1});
    k32 = k32 - L35*(K55{i-1}\R52{i-1});
    k41 = k41 - L45*(K55{i-1}\R51{i-1});
    k42 = k42 - L45*(K55{i-1}\R52{i-1});
    
    b(:,i) = b(:,i) + m(i)*d(:,i)/2;
    c(:,i) = c(:,i) + K33{i-1}\c(:,i-1) - L34*(K44{i-1}\d(:,i-1)) - L35*(K55{i-1}\e(i-1));
    d(:,i) = d(:,i) + K44{i-1}\d(:,i-1) - L45*(K55{i-1}\e(:,i-1));
    
    %upper triangularization of matrix
    k32 = k32 + k31*S(r(:,i));
    k33 = eye(3) + k31*Icm{i};
    k33 = k33 + k32*m(i)*S(l(:,i))/2;
    k34 = k32*m(i);
    k42 = k42 + k41*S(r(:,i));
    k43 = S(l(:,i)) + k41*Icm{i};
    r31 = k31;
    r32 = k31*S(r(:,i));
    r32 = r32 + k32;
    r41 = k41;
    r42 = k41*S(r(:,i));
    k43 = k43 + k42*m(i)*S(l(:,i))/2;
    k44 = eye(3) + k42*m(i);
    r42 = r42 + k42;
    k44 = k44 - k43*(k33\k34);
    k45 = k43*(k33\n(:,i));
    r41 = r41 - k43*(k33\r31);
    r42 = r42 - k43*(k33\r32);
    k53 = n(:,i)'*Icm{i} + n(:,i)'*S(r(:,i))*m(:,i)*S(l(:,i))/2;
    k54 = n(:,i)'*S(r(:,i))*m(i) - k53*(k33\k34);
    k55 = k53*(k33\n(:,i));
    r51 = n(:,i)' - k53*(k33\r31);
    r52 = 2*n(:,i)'*S(r(:,i)) - k53*(k33\r32);
    k55 = k55 - k54*(k44\k45);
    r51 = r51 - k54*(k44\r41);
    r52 = r52 - k54*(k44\r42);
    
    c(:,i) = c(:,i) + k31*a(:,i) + k32*b(:,i);
    d(:,i) = d(:,i) + k41*a(:,i) + k42*b(:,i) - k43*(k33\c(:,i));
    e(i) = e(i) + n(:,i)'*a(:,i) + n(:,i)'*S(r(:,i))*b(:,i) - k53*(k33\c(:,i)) - k54*(k44\d(:,i));
    
    
    K33{i} = k33;
    K34{i} = k34;
    K35{i} = -n(:,i);
    R31{i} = r31;
    R32{i} = r32;
    K44{i} = k44;
    K45{i} = k45;
    R41{i} = r41;
    R42{i} = r42;
    K55{i} = k55;
    R51{i} = r51;
    R52{i} = r52;
    
    
end

%backward recursion
qdd(NB)   = K55{NB}\(e(NB));
acc(:,NB) = K44{NB}\(d(:,NB) - K45{NB}*qdd(NB));
wd(:,NB)  = K33{NB}\(c(:,NB) - K34{NB}*acc(:,NB) + n(:,NB)*qdd(:,NB));
F(:,NB)   = -(b(:,NB) - m(NB)*acc(:,NB) - m(NB)*S(l(:,NB))*wd(:,NB)/2);
T(:,NB)   = -(a(:,NB) - S(r(:,NB))*F(:,NB) - Icm{NB}*wd(:,NB));

for i=NB-1:-1:1
    qdd(i)   = K55{i}\(e(i)                                       - R51{i}*T(:,i+1) - R52{i}*F(:,i+1));
    acc(:,i) = K44{i}\(d(:,i) - K45{i}*qdd(i)                     - R41{i}*T(:,i+1) - R42{i}*F(:,i+1));
    wd(:,i)  = K33{i}\(c(:,i) - K34{i}*acc(:,i) + n(:,i)*qdd(:,i) - R31{i}*T(:,i+1) - R32{i}*F(:,i+1));
    F(:,i)   = -(b(:,i) - m(i)*acc(:,i) - m(i)*S(l(:,i))*wd(:,i)/2                       - F(:,i+1));
    T(:,i)   = -(a(:,i) - S(r(:,i))*F(:,i) - Icm{i}*wd(:,i)            - T(:,i+1)        - S(r(:,i))*F(:,i+1));
end
qdd = qdd';
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
