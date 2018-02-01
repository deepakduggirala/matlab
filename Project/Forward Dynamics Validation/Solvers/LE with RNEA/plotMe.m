x=0:pi/100:pi;
y=sin(x);
y2=cos(x);
hold on
grid on
plot(x,y,'--k','LineWidth',2)
plot(x,y2,':k','LineWidth',1.5)
h_legend = legend('sin(x)','cos(x)','FontSize',50);
set(h_legend,'FontSize',14)
title('sin and cos','FontSize',20)