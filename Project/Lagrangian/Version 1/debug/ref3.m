function ref3(nL)

global r0;
global r;
global m;
global n0;
global tht;
global thtd;
global l0;

p = 1;
elemLen = 0.1*ones(nL,1);
m = p*elemLen;
endTime = 500*10^-3;
timeStep = 10^-3;

tau = [0 0]';
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
n = n0;

l0 = zeros(3,nL);
l0(1,:) = elemLen;           %TODO for any angle


figure()
i=1;
fileId = fopen('bookResults.txt','w');
for t=0:timeStep:endTime
    
    D = getD(nL);
    H = getH(nL);
    C = getC(nL);

    thtdd = inv(D)*(tau-H-C);
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
global tht;
tht2 = tht(2);
D = (10^-3)*[(5/3)+cos(tht2) (1/3)+cos(tht2)/2; (1/3)+cos(tht2)/2  (1/3)];
end
function H = getH(nL)
global tht;
global thtd;
tht2 = tht(2);
thtd2 = thtd(2);
thtd1 = thtd(1);

H = (10^-3)*[-0.5*sin(tht2)*(thtd2)^2 - sin(tht2)*thtd1*thtd2; 0.5*sin(tht2)*(thtd1^2)];
end
function C = getC(nL)
global tht;
tht2 = tht(2);
tht1= tht(1);
C = (10^-2)*9.806*[0.5*cos(tht1) + 0.5*cos(tht1+tht2) + cos(tht1); 0.5*cos(tht1 + tht2)];
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
% Ui = i;
% Uj = j;
% U = R;

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
coords;
plot(coords(1,:),coords(2,:),'o-');
end
function printToFile(fileId,t)
global tht;
fprintf(fileId, 't=%0.5f\ttht = %.5f  %.5f\n',t,tht(1),tht(2));
end