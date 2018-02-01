function qdd = FD_ABA(model)
%except model,q,qd,s all are cell arrays
%unpacking model
s = model.s;
NB = model.NB;
g = model.g;
I = model.I;
q = model.q;
qd = model.qd;
tau = model.tau;

% %pre-allocating space
Xup = cell(1,NB);
v = cell(1,NB);
c = cell(1,NB);
% p = cell(1,NB);
IA = cell(1,NB);
pA = cell(1,NB);
% Ia = cell(1,NB);
% pa = cell(1,NB);
% h = cell(1,NB);
% d = cell(1,NB);
% u = cell(1,NB);
ad = cell(1,NB);
qdd = cell(1,NB);
a = cell(1,NB);

tic
%pass 1
    XJ = rotz(q(1));
    VJ = s*qd(1);
    Xup{1} = XJ*model.Xtree{1};
    v{1} = VJ;
    c{1} = zeros(6,1);
    IA{1} = model.I{1};
    pA{1} = crf(v{1})*(model.I{1}*v{1});
for i=2:NB
    XJ = rotz(q(i));
    VJ = s*qd(i);
    Xup{i} = XJ*model.Xtree{i};
    v{i} = Xup{i}*v{i-1} + VJ;
    c{i} = crm(v{i})*VJ;
    IA{i} = model.I{i};
    pA{i} = crf(v{i})*(model.I{i}*v{i});
end
toc

tic
%pass 2
for i=NB:-1:1
h{i} = IA{i} * s;
d{i} = s'*h{i};
u{i} = tau{i} - s'*pA{i};
if i ~= 1
    Ia = IA{i} - h{i}/d{i}*h{i}';
    pa = pA{i} + Ia*c{i} + h{i} * u{i}/d{i};
    IA{i-1} = IA{i-1} + Xup{i}' * Ia * Xup{i};
    pA{i-1} = pA{i-1} + Xup{i}' * pa;
end
end
toc
tic
%pass 3
for i=1:NB
    if(i==1)
        ad{i} = Xup{i}*(-g) + c{i};
    else
        ad{i} = Xup{i}*a{i-1} + c{i};
    end
    qdd{i} = (u{i} - h{i}'*ad{i})/d{i};
    a{i} = ad{i} + s*qdd{i};
end
toc

qdd = cell2mat(qdd);
end