function plotter(model)

s_ref = 0*10^-3;
t_ref = 10^-3;
e_ref = 1000*10^-3;
t2=s_ref:t_ref:e_ref;

NB = model.NB;

M = dlmread('results.txt','\t');
[k,~] = size(M);
for i=1:k
    q = M(i,2:4);
    r = get_r(model,q);
    for j=1:NB
        CoM_pos(i,:,j) = r(1:3,j)';
    end
end
M2 = dlmread('threelink2D.tab','\t');
for i=1:NB
    ref(:,:,i) = M2(:, 3*i-1:3*i+1);
end
colors = {'r','m','b'};
for i=1:NB
    handle = figure();
    title(strcat('Link ',num2str(i)));
    hold on
    grid on
    for j=1:NB
        plot(M(:,1),CoM_pos(:,j,i),colors{j},'LineWidth',1.5);
        plot(M2(:,1),ref(:,j,i),'c','LineWidth',1.5);
    end
    xlabel('Time (s)')
    ylabel('Position of CoM (m)')
%     legend('x-NE','y-NE','z-NE','Simulink');
%     filename = strcat('Link',num2str(i));
%     filename = strcat(filename,'.jpg');
%     saveas(handle,filename);
end
% for l=1:NB 
%     j=1;
%     max_err = 0;
%     for i=1:k
%         t = t_ref*round(M(i,1)/t_ref);
%         while M2(j,1)<t && j<=(e_ref-s_ref)/t_ref
%             j=j+1;
%         end
%         if j>e_ref/t_ref
%             break;
%         end
%         if M2(j,1) == t
%             err = max(abs(ref(j,:,l) - CoM_pos(i,:,l)));
%             if err>max_err
%                 max_err = err;
%                 t_err = t;
%             end
%         end
%     end
%     fprintf('max abs err of link %d is %f at %f\n',l,max_err,t_err);
% end
end

function r = get_r(model,q)
    r0 = model.r0;
    r0(4,:)=1;
    Xnet = eye(4);
    for i=1:model.NB
        Xnet = Xnet*getR(i,model.n0,q(i),model.l0);
        r(:,i) = Xnet*r0(:,i);
    end
        
end

function T = getR(i,n,tht,l)
    T = zeros(4);
    T(4,4) = 1;
    T(1:3,1:3) = rotMat(n(:,i),tht);
    if i~=1
        T(1:3,4) = l(:,i-1);
    end
    
end

function R = rotMat(n,tht)
n = n/(norm(n));
S = [0,-n(3),n(2);
    n(3),0,-n(1);
    -n(2),n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + (S)*sin(tht);
end