function test()
m = getModelN(50);
N = 1;
tic
for i = 1:N
    qdd = LE_solver(m);
end
elapsedTime = toc;
avgTime = elapsedTime/N
end