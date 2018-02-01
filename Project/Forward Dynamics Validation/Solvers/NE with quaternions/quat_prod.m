function x = quat_prod(a,b)     %Quaternion Conjugate
s1 = a(1);
s2 = b(1);
v1 = a(2:4);
v2 = b(2:4);
x = [s1*s2 - dot(v1,v2); s1*v2 + s2*v1 + cross(v1,v2)];
end