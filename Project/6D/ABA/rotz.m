function X = rotz(q)
c = cos(q);
s = sin(q);
E = [c -s 0
     s c 0
     0 0 1];
X = [E zeros(3)
    zeros(3) E];
end