figure();
grid on
axis equal
axis normal
axis([0 1 0 1])
for tht = 0:0.01:pi/2
    x = cos(tht);
    y = sin(tht);
    
    plot([0 x],[0 y],'o-')
    axis([0 1 0 1])
    pause(0.001)
end
