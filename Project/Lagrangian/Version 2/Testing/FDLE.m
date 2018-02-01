function [qdd,D,H,C] = FDLE(model)
%calculating C(potential energy term)
m_C = model;
m_C.qd = zeros(model.NB,1);
m_C.qdd = zeros(model.NB,1);
C = RNEA(m_C);

%calculating H(velocity term)
m_H = model;
m_H.g=[0 0 0]';
m_H.qdd = zeros(model.NB);
H = RNEA(m_H);

%calculating D(Inertia term)
m_D = model;
m_D.g=[0 0 0]';
m_D.qd = zeros(model.NB,1);
D = zeros(model.NB);
I = eye(model.NB);
for i=1:model.NB
    m_D.qdd = I(:,i);
    D(:,i) = RNEA(m_D);
end
qdd = D\(-H-C);
end