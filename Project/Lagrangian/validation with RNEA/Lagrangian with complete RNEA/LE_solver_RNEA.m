function qdd = LE_solver_RNEA(model_RNEA)
%unpacking model
NB = model_RNEA.NB;
tau = model_RNEA.tau;


model_RNEA.g = [0 0 0]';
H1 = RNEA(model_RNEA)

model_RNEA.qd = [0;0;0];
model_RNEA.g = [0 9.80665 0]';
C1 = RNEA(model_RNEA)

model_RNEA.g = [0 0 0]';
I = eye(NB);
for i=1:NB
    model_RNEA.qdd = I(:,i);
    D1(:,i) = RNEA(model_RNEA);
end
D1
qdd = inv(D1)*(tau-H1-C1);
end