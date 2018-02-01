clear
Sr2 = sym('Sr2');
K13 = sym('K13');
K14 = sym('K14');
K23 = sym('K23');
K24 = sym('K24');
n2 = sym('n2');
n2T = sym('n2T');
R53 = sym('R53');
R54 = sym('R54');
Sr3 = sym('Sr3');
Icm3 = sym('Icm3');
m3 = sym('m3');
A = [1 -2*Sr2 K13 K14 0        0 0 0 0 0;
     0 1      K23 K24 0        0 0 0 0 0;
     0 0       1   0  n2       0 0 -1 0 0;
     0 0       0   1 -(Sr2*n2) 0 0 Sr2+Sr3 -1 0;
     0 0       0   0     1     0 0 R53 R54 0
     -1 0      0   0     0     1 -2*Sr3 Icm3 -m3*Sr3 0
     0 -1      0   0     0     0 1 0 m3 0];

 
 for i = 1:5
     f = -A(6,i);
     A(6,i)=0;
     for j=i+1:10
        A(6,j) = A(6,j) + f*A(i,j);
     end
     A(6,:)
 end
 
 A(6,:)