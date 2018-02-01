function rbi = mcI(m,c,I)
C = skew(c);
rbi = [ I + m*C*C', m*C; m*C', m*eye(3) ];
end