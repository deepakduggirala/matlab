clear
Sr2 = sym('Sr2');
K13 = sym('K13');
K14 = sym('K14');
K23 = sym('K23');
K24 = sym('K24');
n2 = sym('n2');
n2T = sym('n2T');
Sr3 = sym('Sr3');
A = [1 -2*Sr2 K13 K14 0 0 0 0 0 0;
     0 1 K23 K24 0 0 0 0 0 0;
     0 0 1 0 n2 0 0 -1 0 0;
     0 0 0 1 -(Sr2*n2) 0 0 (Sr3+Sr2) -1 0;
     n2T 0 0 0 0 0 0 0 0 0];
 
 for i = 1:4
     f = -A(5,i);
     A(5,i)=0;
     for j=i+1:10
        A(5,j) = A(5,j) + f*A(i,j);
     end
     A(5,:)
 end
 A(5,:)