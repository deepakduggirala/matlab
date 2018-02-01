function [B,b] = gauss(B,b)
[R,C] = size(B);
for row = 2:R
    for col = 1:row-1
        xmult = -B(row,col)/B(col,col);
        B(row,col) = 0;
        for i = col+1:C
            B(row,i) = B(row,i) + xmult*B(col,i);
        end
        b(row) = b(row) + xmult*b(col);
    end
end
end