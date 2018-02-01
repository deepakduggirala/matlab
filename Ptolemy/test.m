function x = test()
select = [1 1 1 1 0 0 0 0 0 0]; %select which body path to plot; 1 = yes, 0 = no
names = {'Sun', 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto'};
x = makeLegend(names,select);
end
function legend_ = makeLegend(names,select)
    legend_ = {};
    j = 1;
    for i=1:10
        if select(i)
            legend_{j} = names(i);
            j = j+1;
        end
    end
end