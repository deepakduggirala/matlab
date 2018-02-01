A = [1	0	0	0	0	0	-1.2254e-09	0	0	0	0	0	0
0	1	0	0	0	0.05	0	-2.0429e-06	0	0	0	0	0
0	0	1	0	-0.05	0	0	0	-2.0429e-06	0	0	0	0
0	0	0	1	0	0	0	0	0	-0.0012254	0	0	0
0	0	0	0	1	0	0	0	0	0	-0.0012254	0	0
0	0	0	0	0	1	0	0	0	0	0	-0.0012254	0
0	0	0	0	0	0	1	0	0	0	0	0	0
0	0	0	0	0	0	0	1	0	0	0	0	0
0	0	0	0	0	0	0	0	1	0	0	0	-1
0	0	0	0	0	0	0	0	0	1	0	0	0
0	0	0	0	0	0	0	0	-0.1	0	1	0	0
0	0	0	0	0	0	0	0.1	0	0	0	1	0
0	0	1	0	0	0	0	0	0	0	0	0	0];

sp_A = sparse(A);
b = [0
     0
     0
     0
-0.0240
     0
     0
     0
     0
     0
     0
     0
     0];
 [N,~] = size(A);
 [rowIndA, colIndA, valA] = find(sp_A);
 j=1;
 for i = 1:size(colIndA)-1
    if colIndA(i) == rowIndA(i)
        y_kk = valA(i);
    end
    if colIndA(i) < rowIndA(i) %lower triangle
        valL(j) = -valA(i)/y_kk;
        rowIndL(j) = rowIndA(i);
        colIndL(j) = colIndA(i);
        j = j+1;
    end
 end