function qdd = solver(q,qd)
l=0.1;
qdd = 3*9.80665*cos(q)/(2*l);
end