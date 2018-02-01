function test()
m = getModelN(100);
N = 10;
tic
for i = 1:N
    qdd = LE_solver(m);
end
elapsedTime = toc;
avgTime = elapsedTime/N
end