function main2(NB,report)

startTime = 0*10^-3;
timeStep = 10^-3;
delT = timeStep;
endTime = 1000*10^-3;

s_ref = 0*10^-3;
t_ref = 10^-3;
e_ref = 1000*10^-3;

refDataFileName = 'threeLink3D.tab';

model = getModel(NB);
i=1;

for t=startTime:timeStep:endTime
    m    = model;
    q1   = m.q;
    qd1  = m.qd;
    qdd1 = LE_solver(m);
    
    q2   = q1  + 0.5*qd1*delT;
    qd2  = qd1 + 0.5*qdd1*delT;
    m.q  = q2;
    m.qd = qd2;
    qdd2 = LE_solver(m);
    
    q3   = q1  + 0.5*qd2*delT;
    qd3  = qd1 + 0.5*qdd2*delT;
    m.q  = q3;
    m.qd = qd3;
    qdd3 = LE_solver(m);
    
    q4   = q1  + qd3*delT;
    qd4  = qd1 + qdd3*delT;
    m.q  = q4;
    m.qd = qd4;
    qdd4 = LE_solver(m);
    
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

for i=1:NB
    hX = figure();
    if(flag==1)
        subplot(3,2,1)
    end
    hold on
    grid on
    plot(t1,CoM_pos(:,1,i),'--b','LineWidth',1.5)
    if(strcmp(refDataFileName,'') == 0)
        plot(t2,ref(:,1,i),'g','LineWidth',1.5)
    end
    title(strcat('x of link',num2str(i)),'FontSize',18)
    h_legend = legend('Lagrangian','Simulink');
    set(h_legend,'FontSize',11)
    xlabel('time (sec)')
    ylabel('Position of CoM (m)')
    
    if(flag==1)
        subplot(3,2,2)
        error = CoM_pos(:,1,i)-ref(:,1,i);
        plot(t1,error,'r','LineWidth',1.5)
        grid on
        title('Error')
        xlabel('time (sec)')
        ylabel('Error (m)')
    end
    
%     hY = figure();
    if(flag==1)
        subplot(3,2,3)
    end
    hold on
    grid on
    plot(t1,CoM_pos(:,2,i),'--b','LineWidth',1.5)
    if(strcmp(refDataFileName,'') == 0)
        plot(t2,ref(:,2,i),'g','LineWidth',1.5)
    end
    title(strcat('y of link',num2str(i)),'FontSize',18)
%     h_legend = legend('Lagrangian','Simulink');
%     set(h_legend,'FontSize',11)
    xlabel('time (sec)')
    ylabel('Position of CoM (m)')
    if(flag==1)
        subplot(3,2,4)
        error = CoM_pos(:,2,i)-ref(:,2,i);
        plot(t1,error,'r','LineWidth',1.5)
        grid on
        title('Error')
        xlabel('time (sec)')
        ylabel('Error (m)')
    end
    
    
%     hZ = figure();
    if(flag==1)
        subplot(3,2,5)
    end
    hold on
    grid on
    plot(t1,CoM_pos(:,3,i),'--b','LineWidth',1.5)
    if(strcmp(refDataFileName,'') == 0)
        plot(t2,ref(:,3,i),'g','LineWidth',1.5)
    end
    title(strcat('z of link',num2str(i)),'FontSize',18)
%     h_legend = legend('Lagrangian','Simulink');
%     set(h_legend,'FontSize',11)
    xlabel('time (sec)')
    ylabel('Position of CoM (m)')
    if(flag==1)
        subplot(3,2,6)
        error = CoM_pos(:,3,i)-ref(:,3,i);
        plot(t1,error,'r','LineWidth',1.5)
        grid on
        title('Error')
        xlabel('time (sec)')
        ylabel('Error (m)')
    end
    
    
    if(report==1)
        filename = strcat('link',num2str(i));
        filename = strcat(filename,'-X');
        filename = strcat(filename,'.jpg');
        saveas(hX,filename);
        
%         filename = strcat('link',num2str(i));
%         filename = strcat(filename,'-Y');
%         filename = strcat(filename,'.jpg');
%         saveas(hY,filename);
%         
%         filename = strcat('link',num2str(i));
%         filename = strcat(filename,'-Z');
%         filename = strcat(filename,'.jpg');
%         saveas(hZ,filename);
    end
end

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