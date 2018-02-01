K14 = {-m1*Sr1};
K
i=1;
q_t = -n(:,i)'*(K14{i} + 2*S(r(:,i))*K24{i});
p_t = -n(:,i)'*(K13{i} + 2*S(r(:,i))*K23{i});
D{i} = (-p_t + q_t*S(r(:,i)))*n(:,i)
if i~=NB
    R53{i} = p_t - q_t*(S(r(:,i)) + S(r(:,i+1)));
    R54{i} = q_t;
end