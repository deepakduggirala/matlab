nL=5;
%Initial orientation of links are in random direction
r  = zeros(nL,3);
for i = 1:nL
    r(i,1) = sign(rand-0.5)*rand;
    r(i,2) = sign(rand-0.5)*rand;
    r(i,3) = sign(rand-0.5)*rand;
end
for i = 1:nL
    r(i,:) = r(i,:)/norm(r(i,:));
end

% thetas are calculated from dotting adjacent links
tht = zeros(nL,1);
tht(1) = 0;
for i=2:nL
    tht(i,:) = dot(r(i-1,:),r(i,:));
end

% axis of rotations are calculated from crossing adjacent links
n = zeros(nL,3);
n(1,:) = [0 0 1];
for i = 2:nL
    n(i,:) = cross(r(i-1,:),r(i,:));
end

% figure();
% hold on;
% axes();
% axis equal;
% 
% coords = zeros(nL+1,3);
% coords(1,:) = [0 0 0];
% for i=2:nL+1
%     coords(i,:) = coords(i-1,:) + r(i-1,:); 
% end
% 
% coords
% plot3(coords(:,1),coords(:,2),coords(:,3),'o-')