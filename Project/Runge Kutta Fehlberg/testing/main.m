function main()
    t_start=0;
    t_end=7;
    step = 10^-3;
    i=1;
    y = zeros((t_end-t_start)/step,1);
    y_prev = 0;
    for t=t_start:step:t_end
        y(i) = RKF5(t,y_prev,step);
        y_prev = y(i);
        i=i+1;
    end
    t=t_start:step:t_end;
    hold on
    grid on
    plot(t,y);
    plot(t,sin(t),'r');
end