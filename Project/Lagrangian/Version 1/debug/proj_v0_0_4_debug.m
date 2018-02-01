function proj_v0_0_4_debug(nL)
global J
global thtd;
global g;
global r;
global n;
global tht;
global dm;
global N0;
global R0;
global thtdd;
global J0;
g= [0 9.8 0];
p = 1;                      %Linear density
elemLen = 0.1*ones(nL,1);   %All element Lengths are taken as 1. Row vector
dm = p*elemLen;             %Mass of links. Row vector
drad = min(elemLen)/100;    %radius of the element
timeStep = 10^-3;           %Time step
endTime = 0.5;

tau = [0]';                 %Torque applied on each link 

N0 = zeros(3,nL);
N0(3,:) = 1;

R0 = zeros(3,nL);
R0(1,:) = -0.5*elemLen;

r = R0;
n = N0;

tht = zeros(nL,1);
thtd = zeros(nL,1);
thtdd = zeros(nL,1);

for i=1:nL
    j = eye(3);
    Iyy = dm(i)*((elemLen(i)^2)/12);
    Ixx = 0;
    Izz = dm(i)*((elemLen(i)^2)/12);
    j(1,1) = (Iyy+Izz-Ixx)/2;
    j(2,2) = (Izz+Ixx-Iyy)/2;
    j(3,3) = (Ixx+Iyy-Izz)/2;
    J0(:,:,i) = j;       %J is a 3d matrix
end

J = J0;

 
fileId = fopen('results.txt','w');

i=1;
figure()
hold on
axis equal
title('-2r')
grid on
xlabel('x-axis')
ylabel('y-axis')
for t=0:timeStep:endTime
    printToFile(fileId,t,nL);
    myPlot(tht,r)
    D = getD(nL);
    H = getH(nL);
    C = getC(nL);
    
    thtdd = inv(D)*(tau - H - C);
    thtdd_store(i) = thtdd;
    thtd = thtd + thtdd*timeStep;
    thtd_store(i) = thtd;
    tht = tht + thtd*timeStep;
    tht_store(i) = tht;
    update_r_n(nL);
    i=i+1;
%     pause()
end
I = 0:timeStep:endTime;

figure()
plot(I,tht_store)
title('tht')
grid on

figure()
plot(I,thtd_store)
title('thtd')
grid on

figure()
plot(I,thtdd_store)
title('thtdd')
grid on
% plot(store_tht(:,1),store_tht(:,2))
end

function D = getD(nL)
global J;
D = zeros(nL,nL);
for i=1:nL
    for k=1:nL
        sum=0;
        for j = max(i,k):nL
            a = getU(j,k);
            b = J(:,:,j);
            c = getU(j,i)';
            sum = sum + trace(a*b*c);
        end
        D(i,k) = sum;
    end
end
end

function H = getH(nL)
global thtd;
H = zeros(nL,1);
for i = 1:nL
    sum = 0;
    for k=1:nL
        for m=1:nL
            sum = sum - get_h(i,k,m,nL)*thtd(k)*thtd(m);
        end
    end
    H(i) = sum;
end
end

function h = get_h(i,k,m,nL)
h=0;
global J;
    for j = max(i,max(k,m)):nL
        a = getU2(j,k,m);
        b = J(:,:,j);
        c = getU(j,i);
        c=c';
        h = h + trace(a*b*c);
    end
end

function C = getC(nL)
% C = zeros(nL,1);
% global g;
% global r;
% global dm;
% for i = 1:nL
%     sum=0;
%     for j=i:nL
%         sum = sum + dm(j)*g*getU(j,i)*r(:,j);
%     end
%     C(i) = sum;
% end
end

function update_r_n(nL)
global r;
global tht;
global N0;
global R0;
global J;
global J0;
    for i = 1:nL
        R = eye(3);
        for j=1:i
            R = R*rotMat(N0(:,i),tht(j));
        end
        r(:,i) = R*R0(:,i);
         J(:,:,i) = R*J0(:,:,i)*R';
    end
end

function R = getU(i,j)  
global tht
global N0;
R = eye(3,3);
for ii=1:i
    if ii == j
        R = R*rotMatD(N0(:,i),tht(ii));
%         x=N0(:,i);
%         y=tht(ii);
    else
        R = R*rotMat(N0(:,i),tht(ii));
    end
end
end

function R = getU2(i,j,k)
global tht
R = [-cos(tht) -sin(tht) 0;
    sin(tht) -cos(tht) 0;
    0 0 0];
end
% function R = getU2(i,j,k)
% global N0
% global tht
% R = eye(3,3);
% for ii=1:i
%     if ii == j
%         R = R*rotMatD(N0(:,i),tht(ii));
%     end
%     
%     if ii == k
%         R = R*rotMatD(N0(:,i),tht(ii));
%     end
%     
%     if ii ~= k && ii ~= j
%         R = R*rotMat(N0(:,i),tht(ii));
%     end
% end
% end

function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,n(3),-n(2);
    -n(3),0,n(1);
    n(2),-n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end

function R = rotMatD(n,tht)
n = n/(norm(n));
S = [0,n(3),-n(2);
    -n(3),0,n(1);
    n(2),-n(1),0];
R = (n)*(n')*(sin(tht)) - eye(3)*(sin(tht)) + (S)*cos(tht);
end

function myPlot(tht,r)
plot(-2*r(1),-2*r(2),'o')
end

function printToFile(fileId,t,nL)
global r;
global n;
global tht;
global thtd;
global thtdd;

for i=1:nL
    fprintf(fileId,'t = %f\tlink = %d\t r = [%f,%f,%f]\t tht=[%f]\t thtd=[%f]\t thtdd=[%f]\n',t,i,r(1,i),r(2,i),r(3,i),tht(i),thtd(i),thtdd(i));
end

end
