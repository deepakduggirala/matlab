function main2(NB,report)

startTime = 0*10^-3;
timeStep = 10^-3;
delT = timeStep;
endTime = 1000*10^-3;

s_ref = 0*10^-3;
t_ref = 10^-3;
e_ref = 1000*10^-3;

refDataFileName = 'threeLink3D.tab';

model = getModel3D(NB);
i=1;

for t=startTime:timeStep:endTime
    m    = model;
    q1   = m.q;
    qd1  = m.qd;
    qdd1 = solver(m)';
    
    q2   = q1  + 0.5*qd1*delT;
    qd2  = qd1 + 0.5*qdd1*delT;
    m.q  = q2;
    m.qd = qd2;
    qdd2 = solver(m)';
    
    q3   = q1  + 0.5*qd2*delT;
    qd3  = qd1 + 0.5*qdd2*delT;
    m.q  = q3;
    m.qd = qd3;
    qdd3 = solver(m)';
    
    q4   = q1  + qd3*delT;
    qd4  = qd1 + qdd3*delT;
    m.q  = q4;
    m.qd = qd4;
    qdd4 = solver(m)';
    
    qf   = q1  + (delT/6)*(qd1 + 2*qd2 + 2*qd3 + qd4);
    qdf  = qd1 + (delT/6)*(qdd1 + 2*qdd2 + 2*qdd3 + qdd4);
    model.q = qf;
    model.qd = qdf;
    r = get_r(model);
    for j=1:NB
        CoM_pos(i,:,j)=r(1:3,j)';
    end
    i=i+1;
end

if(strcmp(refDataFileName,'') == 0)
    M = dlmread(refDataFileName,'\t');
    for i=1:NB
%         ref(:,:,i) = M(:, 3*i-1:3*i+1)./1000;    %1000 is for conversion of mm to m
        ref(:,:,i) = M(:, 3*i-1:3*i+1);
    end
    t2=s_ref:t_ref:e_ref;
    if(size(CoM_pos)==size(ref))
        flag=1;
        error = zeros(NB,1);
        time = zeros(NB,1);
        for i=1:NB
            diff = CoM_pos(:,:,i)-ref(:,:,i);
            std(i,:)=sqrt(var(diff));
            [C,I] = max(abs(diff));
            [C,I1] = max(C);
            error(i) = C;
            time(i) = I(I1); 
        end
        if (report==1)
            fileID = fopen('report.txt','w');
            fprintf(fileID,'No:of links: %d\n',NB);
            fprintf(fileID,'\nsimulation settings:\nendTime = %f\ntimeStep = %f\n',endTime,timeStep);
            fprintf(fileID,'Maximum Absolute Errors:\n');
            for i=1:NB
                fprintf(fileID,'Link-%d\t%f @ %f\n',i,error(i),time(i));
                fprintf(fileID,'Link-%d\t%f\t%f\t%f\n',i,std(i,1),std(i,2),std(i,3));
            end
            fclose(fileID);
        end
    else
        flag=0;
    end
end

t1=startTime:timeStep:endTime;
figure();
hold on;
grid on;
% title('Three link chain in 2D');
% xlabel('Time (s)');
% ylabel('COM (m)');
% % plot(t2,ref(:,1,1),'c','LineWidth',8);
% p11 = plot(t1,CoM_pos(:,1,1),'r--','LineWidth',1.5);
% % plot(t2,ref(:,2,1),'c','LineWidth',8);
% p12 = plot(t1,CoM_pos(:,2,1),'m:','LineWidth',2.5);
% % plot(t2,ref(:,3,1),'c','LineWidth',8);
% p13 = plot(t1,CoM_pos(:,3,1),'b-.','LineWidth',1.5);
% 
% % plot(t2,ref(:,1,2),'y','LineWidth',8);
% p21 = plot(t1,CoM_pos(:,1,2),'r--','LineWidth',1.5);
% % plot(t2,ref(:,2,2),'y','LineWidth',8);
% p22 = plot(t1,CoM_pos(:,2,2),'m:','LineWidth',2.5);
% % plot(t2,ref(:,3,2),'y','LineWidth',4);
% p23 = plot(t1,CoM_pos(:,3,2),'b-.','LineWidth',1.5);
% 
% p31 = plot(t1,CoM_pos(:,1,3),'r--','LineWidth',1.5);
% p32 = plot(t1,CoM_pos(:,2,3),'m:','LineWidth',2.5);
% p33 = plot(t1,CoM_pos(:,3,3),'b-.','LineWidth',1.5);
% 
% color1 = {'r','m',[.49 1 .63]};
% color2 = {'r','m','y'};
% color3 = {'r','m','c'};
% 
% for i=1:3
% tt1 = 1:100:1000;
% val = CoM_pos(tt1,i,1);
% plot(tt1/1000,val,'o','MarkerEdgeColor','k',...
%                       'MarkerFaceColor',color1{i},...
%                       'MarkerSize',14);
% end
% 
% for i=1:3
% tt1 = 1:100:1000;
% val = CoM_pos(tt1,i,2);
% plot(tt1/1000,val,'^','MarkerEdgeColor','k',...
%                       'MarkerFaceColor',color2{i},...
%                       'MarkerSize',10);
% end
% 
% for i=1:3
% tt1 = 1:100:1000;
% val = CoM_pos(tt1,i,3);
% plot(tt1/1000,val,'s','MarkerEdgeColor','k',...
%                       'MarkerFaceColor',color3{i},...
%                       'MarkerSize',10);
% end
% 
% 
% l = legend([p11,p12,p13,p21,p22,p23,p31,p32,p33],{'x1','y1','z1','x2','y2','z2','x3','y3','z3'});
% c = get(l,'Children');
% 
% 
% set(c(25),'Marker','o','MarkerFaceColor',color1{1});
% set(c(22),'Marker','o','MarkerFaceColor',color1{2});
% set(c(19),'Marker','o','MarkerFaceColor',color1{3},'MarkerEdgeColor','k');
% set(c(16),'Marker','^','MarkerFaceColor',color2{1});
% set(c(13),'Marker','^','MarkerFaceColor',color2{2});
% set(c(10),'Marker','^','MarkerFaceColor',color2{3},'MarkerEdgeColor','k');
% set(c(7),'Marker','s','MarkerFaceColor',color1{1});
% set(c(4),'Marker','s','MarkerFaceColor',color1{2});
% set(c(1),'Marker','s','MarkerFaceColor',color1{3},'MarkerEdgeColor','k');

title('Errors','FontSize',20);
ylabel('Error (m)','FontSize',20);
xlabel('Time (s)','FontSize',20);

p11 = plot(t1,CoM_pos(:,1,1)-ref(:,1,1),'r--','LineWidth',1.5);

p12 = plot(t1,CoM_pos(:,2,1)-ref(:,2,1),'m:','LineWidth',2.5);

p13 = plot(t1,CoM_pos(:,3,1)-ref(:,3,1),'b-.','LineWidth',1.5);


p21 = plot(t1,CoM_pos(:,1,2)-ref(:,1,2),'r--','LineWidth',1.5);

p22 = plot(t1,CoM_pos(:,2,2)-ref(:,2,2),'m:','LineWidth',2.5);

p23 = plot(t1,CoM_pos(:,3,2)-ref(:,3,2),'b-.','LineWidth',1.5);

p31 = plot(t1,CoM_pos(:,1,3)-ref(:,1,3),'r--','LineWidth',1.5);
p32 = plot(t1,CoM_pos(:,2,3)-ref(:,2,3),'m:','LineWidth',2.5);
p33 = plot(t1,CoM_pos(:,3,3)-ref(:,3,3),'b-.','LineWidth',1.5);

color1 = {'r','m',[.49 1 .63]};
color2 = {'r','m','y'};
color3 = {'r','m','c'};

for i=1:3
tt1 = 1:100:1000;
val = CoM_pos(tt1,i,1)-ref(tt1,i,1);
plot(tt1/1000,val,'o','MarkerEdgeColor','k',...
                      'MarkerFaceColor',color1{i},...
                      'MarkerSize',14);
end

for i=1:3
tt1 = 1:100:1000;
val = CoM_pos(tt1,i,2)-ref(tt1,i,2);
plot(tt1/1000,val,'^','MarkerEdgeColor','k',...
                      'MarkerFaceColor',color2{i},...
                      'MarkerSize',10);
end

for i=1:3
tt1 = 1:100:1000;
val = CoM_pos(tt1,i,3)-ref(tt1,i,3);
plot(tt1/1000,val,'s','MarkerEdgeColor','k',...
                      'MarkerFaceColor',color3{i},...
                      'MarkerSize',10);
end


l = legend([p11,p12,p13,p21,p22,p23,p31,p32,p33],{'ex1','ey1','ez1','ex2','ey2','ez2','ex3','ey3','ez3'});
% l = legend([p11,p12,p13,p21,p22,p23],{'ex1','ey1','ez1','ex2','ey2','ez2'});
set(l,'FontSize',20);
c = get(l,'Children');


set(c(25),'Marker','o','MarkerFaceColor',color1{1});
set(c(22),'Marker','o','MarkerFaceColor',color1{2});
set(c(19),'Marker','o','MarkerFaceColor',color1{3},'MarkerEdgeColor','k');
set(c(16),'Marker','^','MarkerFaceColor',color2{1});
set(c(13),'Marker','^','MarkerFaceColor',color2{2});
set(c(10),'Marker','^','MarkerFaceColor',color2{3},'MarkerEdgeColor','k');
set(c(7),'Marker','s','MarkerFaceColor',color1{1});
set(c(4),'Marker','s','MarkerFaceColor',color1{2});
set(c(1),'Marker','s','MarkerFaceColor',color1{3},'MarkerEdgeColor','k');

% set(c(16),'Marker','o','MarkerFaceColor',color1{1});
% set(c(13),'Marker','o','MarkerFaceColor',color1{2});
% set(c(10),'Marker','o','MarkerFaceColor',color1{3},'MarkerEdgeColor','k');
% set(c(7),'Marker','^','MarkerFaceColor',color2{1});
% set(c(4),'Marker','^','MarkerFaceColor',color2{2});
% set(c(1),'Marker','^','MarkerFaceColor',color2{3},'MarkerEdgeColor','k');
end




function r = get_r(model)
[R] = preCalc(model.NB,model.l0,model.n0,model.q);

r = zeros(4,model.NB);
model.r0(4,:)=1;
Xnet = eye(4);
for i=1:model.NB
    Xnet = Xnet*R(:,:,i);
    r(:,i) = Xnet*model.r0(:,i);
end
end

function [R] = preCalc(n,l0,n0,q)
R = zeros(4,4,n);
for i = 1:n
    R(:,:,i) = getR(i,l0,n0,q);
end
end

function T = getR(i,l0,n0,tht)
T = zeros(4,4);
T(1:3,1:3) = rotMat(n0(1:3,i),tht(i));
T(4,4) = 1;
if i == 1
    T(1:3,4) = [0 0 0]';
else
    T(1:3,4) = l0(:,i-1);
end
end

function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end