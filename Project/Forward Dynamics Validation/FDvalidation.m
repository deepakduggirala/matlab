function FDvalidation(NB,report)
%Forward Dynamics Validation
%Input:    Forward Dynamics Solver
%          Simulation settings
%          Data to validate (x,y,z) of CoM of each link from solver
%          Data to validate against: ADAMS tab seperated spreadsheet
%          Option to generate report
%                                                  
%
%Features: Works for any model
%          Works for any solver
%          Displays maximum absolute error for each link
%          Displays time of maximum absolute error for each link
%          Three Graphs(x,y,z) for each link along with reference data
%          Generate report on demand which contains graphs as jpgs
%                                                   max abs errors
%
%options to be set at the start
%1.path
%2.reference data filename
%3.simulation time
%4.timeStep
%5.reference data endtime and timestep if decoupled
%
%
%--------------------------------------------------------------------------
%settings
myPath = genpath('C:\Users\deep\Documents\MATLAB\project\Forward Dynamics Validation\Solvers\NE');
refDataFileName = 'L3_2D_CoMpos_e_3_1s.tab';   %filename of reference data

%simulation settings for calculation
endTime = 1000*10^-3;
timeStep = 10^-3;
%reference data timestep and endtime
t_ref = 10^-3;
e_ref = 1;

%DON'T FORGET TO GIVE ROBOT TO getModel

%--------------------------------------------------------------------------
addpath(myPath);
%--------------------------------------------------------------------------
%running solver
numSteps = (endTime/timeStep) +1;
CoM_pos =zeros(numSteps,3,NB);
model = getModel(NB);
i=1;
for t=0:timeStep:endTime
    [qdd,r] = FDNE(model);
    model.qd = model.qd + qdd*timeStep;
    model.q = model.q + model.qd*timeStep;
    for j=1:NB
        CoM_pos(i,:,j)=r(1:3,j)';
    end
    i=i+1;
end
model.q
model.qd

%reading reference data
if(strcmp(refDataFileName,'') == 0)
    M = dlmread(refDataFileName,'\t');
    for i=1:NB
        ref(:,:,i) = M(:,3*i-1:3*i+1)./1000;    %1000 is for conversion of mm to m
    end
    t2=0:t_ref:e_ref;
    if(size(CoM_pos)==size(ref))
        error = zeros(NB,1);
        time = zeros(NB,1);
        for i=1:NB
            diff = CoM_pos(:,:,i)-ref(:,:,i);
            [C,I] = max(abs(diff));
            [C,I1] = max(C);
            error(i) = C
            time(i) = I(I1) 
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

t1=0:timeStep:endTime;

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
rmpath(myPath)
end
