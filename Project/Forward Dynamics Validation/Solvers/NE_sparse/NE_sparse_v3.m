function qdd = NE_sparse_v3(model)
%sparse implimentation of NE_sparse_NE3
% model = getModelSin(100);
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

C = zeros(13*NB,1);
%calculating angular velocities
w(:,1)     = qd(1)*n(:,1);
for i=2:NB
    w(:,i) = qd(i)*n(:,i) + w(:,i-1);
end

%Calculating R.H.S of eqns

for k=1:NB
    i=k;
    a = -S(w(:,i))*(Icm{i}*w(:,i));
    b = m(i)*grav;
    c = S(w(:,i))*(qd(i)*n(:,i));
    d = S(w(:,i))*S(w(:,i))*l(:,i);
    e = tau(i);
    C((k-1)*13+1:k*13,1) = [a;b;c;d;e];
end
C
nnz = 132*NB-54;
I = zeros(nnz,1);
J = zeros(nnz,1);
X = zeros(nnz,1);
kk = 1;
A = {-eye(3), S(r(:,1)), Icm{1}, -eye(3), m(1)*eye(3)/2, eye(3) ,S(l(:,1)), eye(3)};
if NB~=1
    A2 = {eye(3), S(r(:,1)), eye(3)};
    A3 = {m(2)*eye(3)/2, -eye(3), -eye(3)};
end
sI = {[0,0], [0,3], [0,6], [3,3], [3,9], [6,6], [9,6], [9,9]};  %starting Indices
sI2 = {[0,0], [0,3], [3,3]};
sI3 = {[3,9], [6,6], [9,9]};
tic
for k=1:NB
    i=k;
    A{2} = S(r(:,i));
    A{3} = Icm{i};
    A{5} = m(i)*eye(3)/2;
    A{7} = S(l(:,i));
    row = (k-1)*13;
    col = (k-1)*13;
    for ll = 1:8
        x = A{ll};
        temp = sI{ll};
        sub_row = temp(1);
        sub_col = temp(2);
        for ii = 1:3
            for jj=1:3
                I(kk) = row + sub_row + ii;
                J(kk) = col + sub_col +jj;
                X(kk) = x(ii,jj);
                kk = kk+1;
            end
        end
    end
    x = -n(:,i);
    for ii=1:3
        I(kk) = row + 6 + ii;
        J(kk) = col + 13;
        X(kk) = x(ii);
        kk = kk + 1;
    end
    x = n(:,i)';
    for jj=1:3
        I(kk) = row+13;
        J(kk) = col + jj;
        X(kk) = x(jj);
        kk = kk + 1;
    end
    if k<NB
        A2{2} = S(r(:,i));
        row = (k-1)*13;
        col = (k)*13;
        for ll=1:3
            x = A2{ll};
            temp = sI2{ll};
            sub_row = temp(1);
            sub_col = temp(2);
            for ii = 1:3
                for jj=1:3
                    I(kk) = row + sub_row +ii;
                    J(kk) = col + sub_col +jj;
                    X(kk) = x(ii,jj);
                    kk = kk+1;
                end
            end
        end
        A3{1} = m(i+1)*eye(3)/2;
        row = (k)*13;
        col = (k-1)*13;
        for ll=1:3
            x = A3{ll};
            temp = sI3{ll};
            sub_row = temp(1);
            sub_col = temp(2);
            for ii = 1:3
                for jj=1:3
                    I(kk) = row + sub_row +ii;
                    J(kk) = col + sub_col +jj;
                    X(kk) = x(ii,jj);
                    kk = kk+1;
                end
            end
        end
    end
    
end
toc
tic
B = sparse(I,J,X,13*NB,13*NB);
toc
tic
x = B\C;
toc

qdd = zeros(NB,1);
for i=1:NB
    qdd(i) = x(13*i);
end
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
