function U = getU(i,j)
%iInitial axis of rotation is assumed as [0 0 1]'
nL=2;
n = zeros(3,nL);
n(3,:) = 1;
tht = [degtorad(30),degtorad(30)]'
U = eye(3);
for ii = 1:i
    normal = n(:,ii);
    if ii == j
        U = U*rotMatD(normal,tht(ii));
    else
        U = U*rotMat(normal,tht(ii));
    end
end
end

function R = rotMat(n,tht)
%calculates rotational transformation matrix
S = [0,n(3),-n(2);
    -n(3),0,n(1);
    n(2),-n(1),0];
R = (n)*(n')*(1-cos(tht)) + eye(3)*(cos(tht)) + S*sin(tht);
end

function Rd = rotMatD(n,tht)
S = [0,n(3),-n(2);
    -n(3),0,n(1);
    n(2),-n(1),0]';
S2 = [-(n(3)^2)-(n(2)^2) n(1)*n(2) n(1)*n(3)
      n(1)*n(2) -(n(3)^2)-(n(1)^2) n(2)*n(3)
      n(1)*n(3) n(2)*n(3) -(n(2)^2)-(n(1)^2)];
Rd = S2*sin(tht)+S*cos(tht);
end