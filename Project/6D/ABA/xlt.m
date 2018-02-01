function X = xlt(r)
X = [eye(3) zeros(3)
    -skew(r) eye(3)];
end