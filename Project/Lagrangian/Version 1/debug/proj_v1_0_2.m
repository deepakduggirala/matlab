function proj_v1_0_2(nL)
%Working for n=2 and 3 links
%Maximum Error is 0.07883 for n=3 links with time step 10^-4 
global r0;
global r;
global J;
global J0;
global m;
global n0;
global tht;
global thtd;
global l0;

p = 2.4507564291E-004;
elemLen = 0.1*ones(nL,1);
m = p*elemLen;
endTime = 1000*10^-3;
timeStep = 10^-4;

tau = [0 0 0]';
thtd=zeros(nL,1);
tht=zeros(nL,1);
% tht = [degtorad(135)];

r0 = zeros(4,nL);
r0(1,:) = 0.5*elemLen;
r0(4,:) = 1;
r = r0;

n0 = zeros(4,nL);
n0(3,:) = 1;
n0(4,:) = 1;

l0 = zeros(3,nL);
l0(1,:) = elemLen;           %TODO for any angle

for i=1:nL
    j = zeros(4,4);
    Iyy = m(i)*((elemLen(i)^2)/3);
    Ixx = 0;
    Izz = m(i)*((elemLen(i)^2)/3);
    j(1,1) = (Iyy+Izz-Ixx)/2;
    j(2,2) = (Izz+Ixx-Iyy)/2;
    j(3,3) = (Ixx+Iyy-Izz)/2;
    j(:,4) = 0;
    j(4,:) = 0;
    j(4,4) = m(i);
    j(1,4) = m(i)*elemLen(i)/2;
    j(4,1) = m(i)*elemLen(i)/2;
    J0(:,:,i) = j;          %J is a 3d matrix
end
J=J0;

figure()
i=1;
fileId = fopen('myResults.txt','w');
for t=0:timeStep:endTime
    
    D = getD(nL);
    H = getH(nL);
    C = getC(nL);

    thtdd = inv(D)*(tau-H+C);
    thtd = thtd + thtdd*timeStep;
    tht = tht + thtd*timeStep;
    update_r_J(nL);
    if(mod(i,1) ==0)
        myPlot(nL)
    end
    axis equal
    axis normal
    xlabel('X-axis')
    ylabel('Y-axis')
    grid on
    axis([-0.35 0.35 -0.35 0.35])
    pause(0.0001)
    printToFile(fileId,t);
    i=i+1;
end
fclose(fileId);
end

function D=getD(nL)
global J
D = zeros(nL,nL);
for i=1:nL
    for k = 1:nL
        sum = 0;
        for j=max(i,k):nL
            sum = sum + trace(getU(j,k)*J(:,:,j)*(getU(j,i)'));
        end
        D(i,k) = sum;
    end
end
end
function H = getH(nL)
global thtd;
H = zeros(nL,1);
for i=1:nL
    sum = 0;
    for k = 1:nL
        for m = 1:nL
            sum = sum + get_h(nL,i,k,m)*thtd(k)*thtd(m);
        end
    end
    H(i) = sum;
end
end
function sum =get_h(nL,i,k,m)
global J;
sum = 0;
for j=max(i,max(k,m)):nL
    sum = sum + trace(getU2(j,k,m)*J(:,:,j)*(getU(j,i)'));
end
end
function C = getC(nL)
global r0
global m
C = zeros(nL,1);
for i=1:nL
    sum = 0;
    for j=i:nL
      sum = sum + m(j)*[0 -9.8 0 0]*getU(j,i)*r0(:,j);
    end
    C(i) = sum;
end
end
function update_r_J(nL)
global r0;
global r;

for i = 1:nL
    R = eye(4);
    for j=1:i
        R = R*getR(j);
    end
    r(:,i) = R*r0(:,i);
end
end

function R = rotMatDD(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,n(1);
    -n(2),n(1),0];
R = (n*n')*cos(tht) - eye(3)*cos(tht) - S*sin(tht); 
end
function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end
function R = rotMatD(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,n(1);
    -n(2),n(1),0];
R = (n)*(n')*(sin(tht)) - eye(3)*(sin(tht)) + (S)*cos(tht);
end

function T = getR(i)
global l0;
global n0;
global tht;
T = zeros(4,4);
T(1:3,1:3) = rotMat(n0(1:3,i),tht(i));
T(4,4) = 1;
if i == 1
    T(1:3,4) = [0 0 0]';
else
    T(1:3,4) = l0(:,i-1);
end
end
function T = getRD(i)
global n0;
global tht;
T = zeros(4);
T(1:3,1:3) = rotMatD(n0(1:3,i),tht(i));
T(1:3,4) = 0;
T(4,4) = 0;
end
function T = getRDD(i)
global n0;
global tht;
T = zeros(4);
T(1:3,1:3) = rotMatDD(n0(1:3,i),tht(i));
T(1:3,4) = 0;
T(4,4) = 0;
end

function R = getU(i,j)
R = eye(4);
for ii = 1:i
    if ii==j
        R = R*getRD(ii);
    else
        R=R*getR(ii);
    end
end
end
function R = getU2(i,j,k)
R = eye(4);
if i > max(j,k)
    if j==k
        for ii=1:i
            if ii==j
                R = R*getRDD(ii);
            else
                R = R*getR(ii);
            end
        end
    else
        for ii=1:i
            if ii==j
                R = R*getRD(ii);
            end
            if ii==k
                R=R*getRD(ii);
            end
            if ii ~= j && ii ~= k
                R = R*getR(ii);
            end
        end
    end
else
    if i == j && i ~= k
        for ii=1:i
            if ii==j
                R = R*getRDD(ii);
            else
                R = R*getR(ii);
            end
        end
    end
    if i == k && i~= j
        for ii=1:i
            if ii==k
                R = R*getRDD(ii);
            else
                R = R*getR(ii);
            end
        end
    end
    if i == j && i== k
        for ii=1:i
            if ii==j
                R = R*getRDD(ii);
            else
                R = R*getR(ii);
            end
        end
    end
end
end

function myPlot(nL)
global r;
coords = zeros(3,nL+1);
coords(:,1) = [0 0 0]';
for i=2:nL+1
    coords(:,i) = 2*(r(1:3,i-1)) - coords(:,i-1);
end
plot(coords(1,:),coords(2,:),'o-');
end

function printToFile(fileId,t)
global tht;
global r;
fprintf(fileId, '%.10f\n',r(1,1)*1000);
end