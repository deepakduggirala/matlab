function test()
m = getModelN(100);
N = 50;
tic
for i = 1:N
    qdd = NE_rec4(m);
end
elapsedTime = toc;
avgTime = elapsedTime/N
end