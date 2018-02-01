function R = rotMat(n,tht)
%calculates rotational transformation matrix
n = n/(norm(n));
S = [0,n(3),-n(2);
    -n(3),0,n(1);
    n(2),-n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + S*sin(tht);
end