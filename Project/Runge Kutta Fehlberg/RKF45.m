function [k,f] = RKF45(NB)
% m=getModel3D(3); %m is model
% m=getModel(3);
m = getModelSin(NB);
n_max = 1000;   %n_max is maximum number of iterations
t0 = 0; %initial time
tf = 1; %tf is the final time
h0 = 10^-3; %initial timestep

h = h0; 
t = t0; 
k = 0;  %counter for number of steps
f = 0;   %counter for discarded steps

h_min = 10^-5;
h_max = 10;    %min and max bounds for the timestep
e_min = 10^-3;
e_max = 10^-3;    %min and max error tolerance
fileID = fopen('results.txt','w');
tic
while k < n_max && t < tf
    if(h < h_min)
        h = h_min;
    elseif(h > h_max)
        h = h_max;
    end
    [RKF4_q,~,RKF5_q,RKF5_qd] = RKF5(m,h);
    e = max(abs(RKF4_q - RKF5_q));
    if(e>e_max)
        if(h>h_min)
            h = h/2;    %reject the step
            f = f+1;
%             fprintf('step size: %f\n',h);
        end
    else    %accept the step
        k = k+1;
        t = t+h;
        m.q = RKF5_q;
        m.qd = RKF5_qd;
        fprintf(fileID,'%f\t%f\t',t,m.q(1));
        for ii=2:m.NB
            fprintf(fileID,'%f\t',m.q(ii));
        end
        fprintf(fileID,'\n');
        if(e < e_min)
            h = 2*h;
%             fprintf('step size: %f\n',h);
        end
    end
end
toc
fclose(fileID);
% plotter(m);
end