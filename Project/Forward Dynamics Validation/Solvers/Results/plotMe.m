lag = [0.0022,0.1555, 8.6117, 222.1483, 4.4981e+03];
LE_RNEA = [0.0033, 0.0234, 0.1343, 0.4900, 1.8445];
NE_Naive = [0.0015, 0.0047, 0.0118, 0.0234, 0.0465];
NE_featherstone = [0.0014, 0.0042, 0.0099, 0.0192, 0.0380];
NE_sparse = [0.0008, 0.0023, 0.0055, 0.0109, 0.0215];
x = [3,10,25,50,100];
%loglog(x,lag,x,LE_RNEA,x,NE_Naive,x,NE_featherstone,x,NE_sparse);
loglog(x,lag,'-s',x,LE_RNEA,'-+',x,NE_Naive,'-o',x,NE_featherstone,'-d',x,NE_sparse,'-*','LineWidth',1.5);
xlabel('log(No:of links)','FontSize',20);
ylabel('log(Computational time (s))','FontSize',20);
title('Comparision of Forward Dynamics methods','FontSize',20);
set(gca,'FontSize',20);
l = legend('Lagrangian','Lagragian with RNEA', 'NE - Method 1', 'NE - Method 2', 'NE-sparse');
set(l,'FontSize',20);