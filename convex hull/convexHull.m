clear all
hold on
grid on
N=10;
[X,Y] = ginput(N);
plot(X,Y,'o')
for i= 1:N
    points(i).x = X(i);
    points(i).y = Y(i);
end
miny = points(1).y;
minIdx = 1;
for i =2:N
    if(miny > points(i).y)
        miny = points(i).y
        minIdx = i
    elseif(miny == points(i).y)
        if(points(minIdx).x > points(i).x)
            minIdx = i
        end
    end
end
plot(points(minIdx).x, points(minIdx).y, 'ro')
