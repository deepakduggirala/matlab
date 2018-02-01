function vs = crf(v)
    vs = [skew(v(1:3)) zeros(3)
          skew(v(4:6)) skew(v(1:3))];
end