function test()
m = getModelN(100);
N = 1000;
tic
for i = 1:N
    qdd = NE_sparse_v3(m);
end
elapsedTime = toc;
avgTime = elapsedTime/N
end