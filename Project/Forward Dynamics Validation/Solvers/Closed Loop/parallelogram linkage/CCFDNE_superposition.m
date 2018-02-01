function [G,H] = CCFDNE_superposition(model)
NB    = model.NB;
q     = model.q;
qd    = model.qd;
n0    = model.n0;
l0    = model.l0;

G = zeros(6,NB);
H = zeros(6,1);

w = zeros(3,NB);
w(:,1)     = qd(1)*n0(:,1);
for i=2:NB
    w(:,i) = qd(i)*n0(:,i) + w(:,i-1);
end

temp = zeros(3,3);
for i=NB:-1:1
    temp = temp + S(l0(:,i));
    G(1:3,i) = temp*n0(:,i);
    if i>1
        H(1:3,1) = temp*S(w(:,i-1))*qd(i)*n0(:,i);
    end
end

for i=1:NB
    H(1:3,1) = H(1:3,1) + S(w(:,i))*S(w(:,i))*l0(:,i);
end

end