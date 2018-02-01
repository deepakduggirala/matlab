function main()
center = [0, 0];
% Sun = 1;
% Mercury = 2;
% Venus = 3;
% Earth = 4;
% Mars = 5;
% Jupiter = 6;
% Saturn = 7;
% Uranus = 8;
% Neptune = 9;
% Pluto = 10;
radius = [0, 5.79, 10.8, 15, 22.8, 77.8, 143, 287, 450, 590];       %approx semi-major axis
T = [Inf, 0.241, 0.615, 1, 1.88, 11.9, 29.5, 84, 165, 248];
color = ['y','m','c','b','r','g','k',[5,247,25]/256,[191,82,177]/256,[191,144,82]/256];
t = 10;      %time of simulation in years;
N = 10000;   %no:of plot points
select = [1 1 1 1 1 1 1 1 1 1]; %select which body path to plot; 1 = yes, 0 = no
random_pos = 1;     %0 = bodies start on x-axis at start
                    %1 = bodies start on random angles at start

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = 2*pi./T;   %rad/year
if random_pos
    start = 2*pi*rand(1,10);
    start(4) = 0;
else
    start = zeros(1,10);
end

bodies = zeros(3,N,10);
bodies(3,:,:) = 1;
for i=1:10
    [x,y] = circle(radius(i),center,N,w(i),t,start(i));
    bodies(1:2,:,i) = [x;y];
end
plot_c(bodies,color,['Copernian System ', num2str(t), ' years'],select);
plot(0,0,'yo','MarkerFaceColor','y')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bodies_ptolemy = zeros(3,N,10);
bodies_ptolemy(3,:,:) = 1;
Earth = bodies(1:2,:,4);
for i=1:10
bodies_ptolemy(1:2,:,i) = bodies(1:2,:,i) - Earth;
end
plot_c(bodies_ptolemy,color,['Ptolemiac System, ',num2str(t),' years'],select)
plot(0,0,'bo','MarkerFaceColor','b')
end

function plot_c(bodies,color,t,select)
    names = {'Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'};
    figure()
    axis equal
    hold on
    for i = 1:10
        if select(i)
        plot(bodies(1,:,i),bodies(2,:,i),color(i))
        end
    end
    title(t)
    legend(makeLegend(names,select));
end

function [x,y] = circle(r,center,N,w,t,s)
    xc = center(1);
    yc = center(2);
    tht = linspace(s,s+w*t,N);
    x = r*cos(tht) + xc;
    y = r*sin(tht) + yc;
end

function legend_ = makeLegend(names,select)
    legend_ = {};
    j = 1;
    for i=1:10
        if select(i)
        legend_{j} = names{i};
        j = j+1;
        end
    end
end