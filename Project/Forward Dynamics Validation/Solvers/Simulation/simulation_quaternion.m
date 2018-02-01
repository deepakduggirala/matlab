function [Q_solver, Qdd_solver] = simulation_quaternion(model, solver, integrator, refDataFileName, saveIt, planar)
sprintf('Not yet implemented\n');
if strcmp(integrator,'Euler')
    [Q_solver, Qdd_solver] = simulation_quaternion_euler(model, solver,refDataFileName, saveIt, planar);
elseif strcmp(integrator,'RK4')
    sprintf('Not yet implemented\n');
%     simulation_quaternion_RK4(model, solver, report, refDataFileName, saveIt, planar);
else
    sprintf('Unknown solver specified\n');
end
end

function [Q_solver,Qdd_solver] = simulation_quaternion_euler(model_arg, solver, refDataFileName, saveIt, planar)
%simulation settings for calculation
endTime = 1000*10^-3;
timeStep = 10^-3;

%copying model
model = model_arg;
%running solver
numSteps = (endTime/timeStep) +1;
CoM_pos =zeros(numSteps,3,model.NB);
i=1;
%start - debug
Q_solver = zeros(numSteps,1+4*model.NB);
Qdd_solver = zeros(numSteps, 1+4*model.NB);
%end - debug
for t=0:timeStep:endTime
    [Qdd,r] = solver(model);
    Qdf = model.Qd + timeStep*Qdd;
    Qf = zeros(4,model.NB);
    for k=1:model.NB
        Qf(:,k) = getUpdatedQuat(Qdf(:,k), model.Q(:,k), timeStep);
    end
    %start - debug
    Q_flat = zeros(1,4*model.NB);
    for k=1:model.NB
        start_index = 1+(k-1)*4;
        end_index = 4+(k-1)*4;
        Q_flat(start_index:end_index) = Qf(:,k)';
    end
    Q_solver(i,:) = [t,Q_flat];
    Qdd_flat = zeros(1,4*model.NB);
    for k=1:model.NB
        start_index = 1+(k-1)*4;
        end_index = 4+(k-1)*4;
        Qdd_flat(start_index:end_index) = Qdd(:,k)';
    end
    Qdd_solver(i,:) = [t,Qdd_flat];
    %end - debug
    %updating model
    model.Q = Qf;
    model.Qd = Qdf;
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
    
    if(planar == 0)
    hZ = figure();
    hold on
    grid on
    plot(t1,CoM_pos(:,3,i),'r')
    if(strcmp(refDataFileName,'') == 0)
        plot(t2,ref(:,3,i),'b')
    end
    title(strcat('z of link',num2str(i)))
    legend(strcat('z-',func2str(solver)),'z-Simulink')
    end

    if(saveIt==1)
        filename = strcat('link',num2str(i));
        filename = strcat(filename,'-X');
        filename = strcat(filename,'.jpg');
        saveas(hX,filename);
        
        filename = strcat('link',num2str(i));
        filename = strcat(filename,'-Y');
        filename = strcat(filename,'.jpg');
        saveas(hY,filename);
        
        if(planar == 0)
        filename = strcat('link',num2str(i));
        filename = strcat(filename,'-Z');
        filename = strcat(filename,'.jpg');
        saveas(hZ,filename);
        end
    end
end

end

function Q_new = getUpdatedQuat(Qd, Q, del_t)
    W = 2*quat_prod(Qd,qConj(Q));
    w = W(2:4);
    w_n = norm(w);
    del_Q = [cos(w_n*del_t/2);sin(w_n*del_t/2)*w/w_n];
    Q_new = quat_prod(del_Q, Q);
end

function x = quat_prod(a,b)     %Quaternion Conjugate
s1 = a(1);
s2 = b(1);
v1 = a(2:4);
v2 = b(2:4);
x = [s1*s2 - dot(v1,v2); s1*v2 + s2*v1 + cross(v1,v2)];
end

function Q1 = qConj(Q)          %Quaternion conjugate
Q1 = [Q(1);-Q(2);-Q(3);-Q(4)];
end

function r_global = get_r(model,r)
r_global = zeros(3,model.NB);
r_global(:,1) = r(:,1);
for i=2:model.NB
    r_global(:,i) = r_global(:,i-1) + r(:,i-1) + r(:,i);
end
end