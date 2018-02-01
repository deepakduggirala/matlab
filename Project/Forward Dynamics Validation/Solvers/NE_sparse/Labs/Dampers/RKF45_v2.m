function [k,f] = RKF45_v2(n)
m=getModelN(n);    %m is model
tol = 10^-3;    %error tolerance
n_max = 10000;   %max number of iterations
t0 = 0;         %initial time
tf = 1;         %tf is the final time
h = 10^-3;     %initial timestep

t = t0; 
k = 0;          %counter for number of steps
f = 0;          %counter for discarded steps

fileID = fopen('results.txt','w');
tic
while k < n_max && t < tf
    [RKF4_q,~,RKF5_q,RKF5_qd] = RKF5(m,h);
    err = max(abs(RKF4_q - RKF5_q)); %taking the maximum of all errors may not be correct
    if err == 0
        s=1;
    else
       s = ((tol)/(2*err))^0.25;
    end
    
    if(err > tol)
        f = f+1;
    else
        k = k + 1;
        t = t + h;
        m.tau(1) = 1*sin(2*pi*t);
        m.q = RKF5_q;
        m.qd = RKF5_qd;
        fprintf(fileID,'%f\t%f\t',t,m.q(1));
        for ii=2:m.NB
            fprintf(fileID,'%f\t',m.q(ii));
        end
        fprintf(fileID,'\n');
    end
%     s
    h = s*h;
end
toc
fclose(fileID);
%plotter(m);
end