function [k,f] = simulator(n)
    tic
    [k,f] = RKF45_v2(n);
    toc
%     M = dlmread('results.txt');
%     q = M(:,2:end);
%     x = size(q);
%     x = x(1);
%     m = getModelN(n);
%     for i=1:x
%         m.q = q(i,:);
%         l = get_l(m);
%         plotChain(l,n);
%         axis equal
%         axis normal
%         xlabel('X-axis','FontSize',20)
%         ylabel('Y-axis','FontSize',20)
%         grid on
%         axis([-1 11 -5 5])
% %         Movie(i)=getframe(gcf);
%         pause(0.01)
%     end
%     movie2avi(Movie,'CahinMovie.avi');
end

function plotChain(l,N)
    coords = zeros(3,N+1);
    coords(:,1) = zeros(3,1);
    coords(:,2:end) = l(1:3,:);
    set(gca,'FontSize',20);
    plot(coords(1,:),coords(2,:),'LineWidth',1.5);
end

function l = get_l(model)
[R] = preCalc(model.NB,model.l0,model.n0,model.q);

l = zeros(4,model.NB);
model.l0(4,:)=1;
Xnet = eye(4);
for i=1:model.NB
    Xnet = Xnet*R(:,:,i);
    l(:,i) = Xnet*model.l0(:,i);
end
end

function [R] = preCalc(n,l0,n0,q)
R = zeros(4,4,n);
for i = 1:n
    R(:,:,i) = getR(i,l0,n0,q);
end
end

function T = getR(i,l0,n0,tht)
T = zeros(4,4);
T(1:3,1:3) = rotMat(n0(1:3,i),tht(i));
T(4,4) = 1;
if i == 1
    T(1:3,4) = [0 0 0]';
else
    T(1:3,4) = l0(:,i-1);
end
end

function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end