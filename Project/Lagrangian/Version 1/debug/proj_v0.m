function proj_v0(nL)
%
%nL is number of links

global r;       % unit vectors along links
global n;       % unit vector - rotation axis
global tht;     % theta - angle of rotation about rotaional axis
global thtd;    % theta dot - angular speed about rotaional axis
global thtdd;   % theta double dot - angular acceleration about rotaional axis
global dm;      % mass of a link
global J;       % moment of inerta of a link in native coordinate system
global delT;    %time step
%radius of the element(cylinder)
drad = 0.01;
%length of the element(cylinder)
elemLength = 1;
%mass of an element
dm = 0.001;
%time step
delT = 0.01;

%Initial orientation of links is assumed as [0 1 0]
r = zeros(nL,3);
r(:,2) = 1;

%iInitial axis of rotation is assumed as [0 0 1]
n = zeros(nL,3);
n(:,3) = 1;

%all rotations are assumed to be 0 radians
tht = zeros(nL,1);

%all velocities are assumed to be 0
thtd = zeros(nL,1);

%generating an empty vector to store acceleration
thtdd = zeros(nL,1);

%calculating moment of inertia for an element
J = eye(3);
Ixx = 0.25*dm*(drad^2) + (dm*(elemLength^2))/3;
Iyy = 0.5*dm*(drad^2);
Izz = 0.25*dm*(drad^2) + (dm*(elemLength^2))/3;
J(1,1) = (Iyy+Izz-Ixx)/2;
J(2,2) = (Izz+Ixx-Iyy)/2;
J(3,3) = (Ixx+Iyy-Izz)/2;

fileId = fopen('results.txt','w');

for t = 0:delT:0.01
    for i = 1:nL
        Dsolve(getTau(t,i,nL),i,nL);
    end
    update_r_and_n(nL);
    printToFile(fileId,t,nL);
end
fclose(fileId);
end

function Dsolve(tau,i,nL)
global tht;
global thtd;
global thtdd;
global delT;
D = getD(i,nL);
C = getC(i,nL);
V = getV(i,nL);
thtdd(i) = (tau-C-V)/D;
thtd(i) = thtd(i) + thtdd(i)*delT;
tht(i)  = tht(i)  + thtd(i)* delT;
end

function D = getD(i,nL)
global J
D = 0;
for j = i:nL
    for k = 1:j
        u1 = getU(j,k);
        u2 = getU(j,i);
        u2 = u2';
        D = D+trace(u1*J*u2);
    end
end
end

function C = getC(i,nL)
global thtd;
global J;
C = 0;
for j=i:nL
    for k=1:j
        for m=1:j
            u1 = getU(j,k);
            u2 = getU(j,i);
            u2 = u2';
            thtd1 = thtd(k,:);
            thtd1 = thtd1';
            thtd2 = thtd(m,:);
            thtd2 = thtd2';
            
            C = C+trace(u1*J*u2*thtd1*thtd2);
        end
    end
end
end

function V = getV(i,nL)
global dm;
global r;
V = 0;
for j = i:nL
    temp = r(j,:);
    temp = temp';
    V = V + getU(j,i)*temp;
end
k = [0 0 1];
V = -dm*V;
V = k*V;
end

function R = rotMat(n,tht)
%calculates rotational transformation matrix

S = [0,n(3),-n(2);
    -n(3),0,n(1);
    n(2),-n(1),0];
R = n*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + S*sin(tht);
end

function Rd = rotMatD(n,tht)
S = [0,n(3),-n(2);
    -n(3),0,n(1);
    n(2),-n(1),0];
Rd = n*(n')*(sin(tht)) - eye(3)*(sin(tht)) + S*cos(tht);
end

function U = getU(i,j)
global n;
global tht;
U = eye(3);
for ii = 1:i
    normal = n(ii,:);
    normal = normal';
    if ii == j
        U = U*rotMatD(normal,tht(ii));
    else
        U = U*rotMat(normal,tht(ii));
    end
end
end

function tau = getTau(t,i, nL)
%torque as a function of time and on the last link
tau = 0;
if i == nL
    tau = 2*t;
end
end

function update_r_and_n(nL)
global tht;
global r;
global n;
for i = 1:nL
    R = eye(3);
    for j = 1:i
        R = R*rotMat(n(i,:),tht(i));
    end
    r1 = r(i,:);
    r1 = r1';
    n1 = n(i,:);
    n1 = n1';
    r(i,:) = (R*r1)';
    r(i,:) = r(i,:)/norm(r(i,:));
    n(i,:) = (R*n1)';
    n(i,:) = n(i,:)/norm(n(i,:));
end
end

function printToFile(fileId,t,nL)
global r;
global n;
global tht;

for i=1:nL
    fprintf(fileId,'t = %4.2f\tlink = %d\t r = [%4.2f,%4.2f,%4.2f]\t n=[%4.2f,%4.2f,%4.2f]\t tht=[%4.2f]\n',t,i,r(i,1),r(i,2),r(i,3),n(i,1),n(i,2),n(i,3),tht(i));
end

end