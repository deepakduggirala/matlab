tht = deg2rad(60);
l = 0.1;
g = -9.81;
c = cos(tht);
s = sin(tht);
sigma = (40*c^2 + 68 + 6*s^4 - 24*(s^2)*(c^2) + 6*(s^3)*(c^2) - s^2)/(24*c);
alpha1 = g/(l*sigma);
alpha2 = g*cos(tht)/(l*sigma);
alpha3 = -2*g*(s^2)/(l*sigma);