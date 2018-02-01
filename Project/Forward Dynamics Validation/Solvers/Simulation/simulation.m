function simulation(model, solver, integrator, refDataFileName, saveIt)
if strcmp(integrator,'Euler')
    simulation_euler(model, solver,refDataFileName, saveIt);
elseif strcmp(integrator,'RK4')
    simulation_RK4(model, solver, report, refDataFileName, saveIt);
else
    sprintf('Unknown solver specified\n');
end
end

function simulation_euler(model, solver, refDataFileName, saveIt)
%simulation settings for calculation
endTime = 1000*10^-3;
timeStep = 10^-3;

%running solver
numSteps = (endTime/timeStep) +1;
CoM_pos =zeros(numSteps,3,model.NB);
i=1;
for t=0:timeStep:endTime
    [qdd,r] = solver(model);
    qdf = model.qd + timeStep*qdd;
    qf = model.q + timeStep*qdf;
    
    %updating model
    model.q = qf;
    model.qd = qdf;
    CoM_pos(i,:,:) = get_r(model,r);
    i=i+1;
end
t1=0:timeStep:endTime;

%reding reference data
if(strcmp(refDataFileName,'') == 0)     %if refDataFileName is not empty
    M = dlmread(refDataFileName,'\t');
    s = size(M);
    ref = zeros(s(1),3,model.NB);     %reference data matrix
    for i=1:model.NB
        ref(:,:,i) = M(:,3*i-1:3*i+1);
    end
    t2=M(:,1);
end

for i=1:model.NB
    hX = figure();
    hold on
    grid on
    plot(t1,CoM_pos(:,1,i),'r')
    if(strcmp(refDataFileName,'') == 0)
        plot(t2,ref(:,1,i),'b')
    end
    title(strcat('x of link',num2str(i)))
    legend(strcat('x-',func2str(solver)),'x-Simulink')
    
    hY = figure();
    hold on
    grid on
    plot(t1,CoM_pos(:,2,i),'r')
    if(strcmp(refDataFileName,'') == 0)
        plot(t2,ref(:,2,i),'b')
    end
    title(strcat('y of link',num2str(i)))
    legend(strcat('y-',func2str(solver)),'y-Simulink')
    
%     hZ = figure();
%     hold on
%     grid on
%     plot(t1,CoM_pos(:,3,i),'r')
%     if(strcmp(refDataFileName,'') == 0)
%         plot(t2,ref(:,3,i),'b')
%     end
%     title(strcat('z of link',num2str(i)))
%     legend(strcat('z-',func2str(solver)),'z-Simulink')

    if(saveIt==1)
        filename = strcat('link',num2str(i));
        filename = strcat(filename,'-X');
        filename = strcat(filename,'.jpg');
        saveas(hX,filename);
        
        filename = strcat('link',num2str(i));
        filename = strcat(filename,'-Y');
        filename = strcat(filename,'.jpg');
        saveas(hY,filename);
        
%         filename = strcat('link',num2str(i));
%         filename = strcat(filename,'-Z');
%         filename = strcat(filename,'.jpg');
%         saveas(hZ,filename);
    end
end
end

function r_global = get_r(model,r)
r_global = zeros(3,model.NB);
r_global(:,1) = r(:,1);
for i=2:model.NB
    r_global(:,i) = r_global(:,i-1) + r(:,i-1) + r(:,i);
end
end
