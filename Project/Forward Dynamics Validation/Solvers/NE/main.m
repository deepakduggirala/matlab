function main(NB,report)

startTime = 0*10^-3;
timeStep = 10^-3;
endTime = 1000*10^-3;

s_ref = 0*10^-3;
t_ref = 10^-3;
e_ref = 1000*10^-3;

refDataFileName = 'twoLink2D_yzplane.tab';

model = getModel(NB);
i=1;

for t=startTime:timeStep:endTime
    [qdd] = FDNE2(model);
    if(isnan(qdd))
        return
    end
    qd_start = model.qd;
    model.qd = model.qd + qdd*timeStep;
    model.q = model.q + qd_start*timeStep + qdd*(timeStep^2)/2;
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
        error = zeros(NB,1);
        time = zeros(NB,1);
        for i=1:NB
            diff = CoM_pos(:,:,i)-ref(:,:,i);
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
            end
            fclose(fileID);
        end
    end
end

t1=startTime:timeStep:endTime;

for i=1:NB
    hX = figure();
    hold on
    grid on
    plot(t1,CoM_pos(:,1,i),'r')
    if(strcmp(refDataFileName,'') == 0)
        plot(t2,ref(:,1,i),'c')
    end
    title(strcat('x of link',num2str(i)))
    legend('x-lag','x-adams')
    
    hY = figure();
    hold on
    grid on
    plot(t1,CoM_pos(:,2,i),'r')
    if(strcmp(refDataFileName,'') == 0)
        plot(t2,ref(:,2,i),'c')
    end
    title(strcat('y of link',num2str(i)))
    legend('y-lag','y-adams')
    
    hZ = figure();
    hold on
    grid on
    plot(t1,CoM_pos(:,3,i),'r')
    if(strcmp(refDataFileName,'') == 0)
        plot(t2,ref(:,3,i),'c')
    end
    title(strcat('z of link',num2str(i)))
    legend('z-lag','z-adams')
    
    if(report==1)
        filename = strcat('link',num2str(i));
        filename = strcat(filename,'-X');
        filename = strcat(filename,'.jpg');
        saveas(hX,filename);
        
        filename = strcat('link',num2str(i));
        filename = strcat(filename,'-Y');
        filename = strcat(filename,'.jpg');
        saveas(hY,filename);
        
        filename = strcat('link',num2str(i));
        filename = strcat(filename,'-Z');
        filename = strcat(filename,'.jpg');
        saveas(hZ,filename);
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