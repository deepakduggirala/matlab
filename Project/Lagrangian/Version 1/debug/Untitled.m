s1=sym('s1');
c1 = sym('c1');

s2 = sym('s2');
c2 = sym('c2');

R2 = [c2 -s2 0
     s2 c2 0
       0 0 1];
R1 = [c1 -s1 0
       s1 c1 0
       0 0 1];
   dR1 = [-s1 -c1 0
            c1 -s1 0
            0 0 0];
 dR2 = [-s2 -c2 0
            c2 -s2 0
            0 0 0];
U11 = dR1;
U21 = dR1*R2;
U22 = R1*dR2;

U11t = [-s1 c1 0;-c1 s1 0;0 0 0];
U22t = [- c1*s2 - c2*s1,c1*c2 - s1*s2,0;s1*s2 - c1*c2, - c1*s2 - c2*s1, 0;0 0 0];
U21t = [- c1*s2 - c2*s1, c1*c2 - s1*s2, 0; s1*s2 - c1*c2, - c1*s2 - c2*s1, 0; 0 0 0];

k = sym('k');
J = [k 0 0;0 0 0;0 0 0];

D11 = trace(U11*J*U11t) + trace(U21*J*(U21t));
D12 = trace(U22*J*U21t);
D21 = trace(U21*J*U22t);
D22 = trace(U22*J*U22t);

X = [  1/k, -1/k;-1/k,  2/k];

