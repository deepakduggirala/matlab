function proj_v0_0_42_debug(nL)
global r0;
global r;
global J;
global J0;
global m;
global n0;
global tht;
global thtd;

p=1;
elemLen = 1*ones(nL,1);
m = p*elemLen;
endTime = 500*10^-3;
timeStep = 10^-4;

tau = zeros(nL,1);
thtd=zeros(nL,1);
% tht=zeros(nL,1);
tht = [degtorad(0),degtorad(45)]';

r0 = zeros(3,nL);
r0(1,:) = 0.5*elemLen;
r = r0;

n0 = zeros(3,nL);
n0(3,:) = 1;
n = n0;

for i=1:nL
    j = eye(3);
    Iyy = m(i)*((elemLen(i)^2)/12);
    Ixx = 0;
    Izz = m(i)*((elemLen(i)^2)/12);
    j(1,1) = (Iyy+Izz-Ixx)/2;
    j(2,2) = (Izz+Ixx-Iyy)/2;
    j(3,3) = (Ixx+Iyy-Izz)/2;
    J0(:,:,i) = j;       %J is a 3d matrix
end
J=J0;

figure()
axis equal
grid on
hold on
xlabel('X-axis')
ylabel('Y-axis')
i=1;
for t=0:timeStep:endTime
    if mod(i,100) == 0
        myPlot(nL)
    end
%     myPlot(nL);
    D = getD(nL);
    H = getH(nL);
    C = getC(nL);
    x = inv(D);
    thtdd = x*(tau-H+C);
    thtd = thtd + thtdd*timeStep;
    tht = tht + thtd*timeStep;
    update_r_J(nL);
    i=i+1;
end
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
      sum = sum + m(j)*[0 -9.8 0]*getU(j,i)*r0(:,j);
    end
    C(i) = sum;
end
end
function update_r_J(nL)
global r0;
global r;
global J;
global J0;
global n0;
global tht;
for i = 1:nL
    R = eye(3);
    for j=1:i
        R = R*rotMat(n0(:,j),tht(j));
    end
    r(:,i) = R*r0(:,i);
    J(:,:,i) =R*J0(:,:,i)*(R');
end
end
function R = getU(i,j)
global n0;
global tht;
R = eye(3);
for ii = 1:i
    if ii==j
        R = R*rotMatD(n0(:,ii),tht(ii));
    else
        R=R*rotMat(n0(:,ii),tht(ii));
    end
end
end
function R = getU2(i,j,k)
global n0;
global tht;
R = eye(3);
if i > max(j,k)
    if j==k
        for ii=1:i
            if ii==j
                R = R*rotMatDD(n0(:,ii),tht(ii));
            else
                R = R*rotMat(n0(:,ii),tht(ii));
            end
        end
    else
        for ii=1:i
            if ii==j
                R = R*rotMatD(n0(:,ii),tht(ii));
            end
            if ii==k
                R=R*rotMatD(n0(:,ii),tht(ii));
            end
            if ii ~= j && ii ~= k
                R = R*rotMat(n0(:,ii),tht(ii));
            end
        end
    end
else
    if i == j && i ~= k
        for ii=1:i
            if ii==j
                R = R*rotMatDD(n0(:,ii),tht(ii));
            else
                R = R*rotMat(n0(:,ii),tht(ii));
            end
        end
    end
    if i == k && i~= j
        for ii=1:i
            if ii==k
                R = R*rotMatDD(n0(:,ii),tht(ii));
            else
                R = R*rotMat(n0(:,ii),tht(ii));
            end
        end
    end
    if i == j && i== k
        for ii=1:i
            if ii==j
                R = R*rotMatDD(n0(:,ii),tht(ii));
            else
                R = R*rotMat(n0(:,ii),tht(ii));
            end
        end
    end
end
end
function R = rotMatDD(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
R = (n*n')*cos(tht) - eye(3)*cos(tht) - S*sin(tht);
end
function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end
function R = rotMatD(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
R = (n)*(n')*(sin(tht)) - eye(3)*(sin(tht)) + (S)*cos(tht);
end
function myPlot(nL)
global r;
coords = zeros(3,nL+1);
coords(:,1) = [0 0 0]';
for i = 2:nL+1
    coords(:,i) = coords(:,i-1) + 2*r(:,i-1);
end
plot(coords(1,:),coords(2,:),'o-')

end
