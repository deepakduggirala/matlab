function proj_v0_0_3(nL)
global tht;
global r;
global n;
global J;
global thtd;
global dm;
%Configuration
p = 1;                  %Linear density
elemLen = ones(nL,1);   %All element Lengths are taken as 1. Row vector
dm = p*elemLen;         %Mass of links. Row vector
drad = 0.01;            %radius of the element
delT = 0.01;            %Time step
tau = [0 0]';           %Torque applied on each link 


%Initial orientation of links is assumed as [-elemlen/2 0 0]'
r = zeros(3,nL);
for i=1:nL
    r(:,i) = [-0.5*elemLen(i);0;0];
end

%iInitial axis of rotation is assumed as [0 0 1]'
n = zeros(3,nL);
n(3,:) = 1;

tht = [degtorad(30),degtorad(30)]'

%all velocities are assumed to be 0
thtd = zeros(nL,1);


%calculating Moment of inertia
for i=1:nL
    j = eye(3);
    Iyy = 0.25*dm(i)*(drad^2) + (dm(i)*(elemLen(i)^2))/3;
    Ixx = 0.5*dm(i)*(drad^2);
    Izz = 0.25*dm(i)*(drad^2) + (dm(i)*(elemLen(i)^2))/3;
    j(1,1) = (Iyy+Izz-Ixx)/2;
    j(2,2) = (Izz+Ixx-Iyy)/2;
    j(3,3) = (Ixx+Iyy-Izz)/2;
    J(:,:,i) = j;         %J is a 3d matrix
end
J
% fileId = fopen('results.txt','w');
% figure();
% axis equal;
% hold on;

for t = 0:delT:0.01
%     printToFile(fileId,t,nL);
%     myPlot(nL, t)
%     pause(0.1)
    D = getD(nL)
    H = getH(nL);
    C = getC(nL);
    right = (tau-H-C);
    left = inv(D);
    thtdd = left*right;
    thtd = thtd + thtdd*delT;
    tht = mod((tht + thtd*delT),pi);
    update_r_n(nL);
    
    %myPlot(nL);
end
% fclose(fileId);
end

function D = getD(nL)
global J;
D = zeros(nL,nL);
for i = 1:nL
    for k = 1:nL
        sum = 0;
        for j = max(i,k):nL
            j
            k
            a = getU(j,k)
            b = J(:,:,j)
            j
            i
            c = getU(j,i)
            x = a*b*((c)')
            sum = sum + trace(x)         %can U be a premultiplied matrix?
        end
        i
        k
        D(i,k) = sum
    end
end
end

function H = getH(nL)
global thtd;
H = zeros(nL,1);
for i=1:nL
    for k=1:nL
        for m=1:nL
            h = getA(i,k,m,nL)*thtd(k)*thtd(m);    %A is h_ikm
        end
    end
    H(i) = h;
end
end

function sum = getA(i,k,m,nL)
global J;
sum = 0;
for j = max(i,max(k,m)):nL
    sum = sum+trace(getU2(j,k,m)*J(:,:,j)*(getU(j,i))');
end
end

function U = getU2(j,k,m)
global n;
global tht;
U = eye(3);
for ii = 1:j
    normal = n(:,ii);
    if ii == k || m
        U = U*rotMatD(normal,tht(ii));
    else
        U = U*rotMat(normal,tht(ii));
    end
end
end

function C = getC(nL)
global r;
global dm;
g = [0 -9.81 0];
C = zeros(nL,1);
for i = 1:nL
    sum=0;
    for j = i:nL
        sum = sum + dm(i)*getU(j,i)*r(:,j);
    end
    C(i) = g*sum;
end  
end

function update_r_n(nL)
global tht;
global r;
global n;
for i = 1:nL
    R = eye(3);
    for j = 1:i
        R = R*rotMat(n(:,i),tht(i));
    end
    r(:,i) = (R*r(:,i));
    r(:,i) = r(:,i)/norm(r(:,i));
    n(:,i) = (R*n(:,i));
    n(:,i) = n(:,i)/norm(n(:,i));
end
end

function R = rotMat(n,tht)
%calculates rotational transformation matrix
S = [0,n(3),-n(2);
    -n(3),0,n(1);
    n(2),-n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + S*sin(tht);
end

function Rd = rotMatD(n,tht)
S = [0,n(3),-n(2);
    -n(3),0,n(1);
    n(2),-n(1),0]';
Rd = n*(n')*(sin(tht)) - eye(3)*(sin(tht)) + S*cos(tht);
end

function U = getU(i,j)
global n;
global tht;
U = eye(3);
for ii = 1:i
    normal = n(:,ii);
    if ii == j
        U = U*rotMatD(normal,tht(ii));
    else
        U = U*rotMat(normal,tht(ii));
    end
end
end

function printToFile(fileId,t,nL)
global r;
global n;
global tht;

for i=1:nL
    fprintf(fileId,'t = %4.2f\tlink = %d\t r = [%4.2f,%4.2f,%4.2f]\t n=[%4.2f,%4.2f,%4.2f]\t tht=[%4.2f]\n',t,i,r(1,i),r(2,i),r(3,i),n(1,i),n(2,i),n(3,i),tht(i));
end

end

function myPlot(nL, t)
global r;
% figure();
axis equal;
hold on;

coords = zeros(3,nL+1);
coords(:,1) = [0 0 0];
for i=2:nL+1
    coords(:,i) = coords(:,i-1) + r(:,i-1); 
end
plot(coords(1,:),coords(2,:),'o-')
% for i=1:nL+1
%     str = num2str(i);
%     text(coords(1,i),coords(2,i),str)
% end
end
