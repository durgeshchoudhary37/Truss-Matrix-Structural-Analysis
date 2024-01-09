% fuction for element stiffness
function K=ElementStiffness( x1,y1,z1,x2,y2,z2,A ,E)
l=sqrt((x2-x1)^2+(y2-y1)^2 + (z2-z1)^2);
a=(x2-x1)/l;
b=(y2-y1)/l;
c=(z2-z1)/l;
T=[a b c 0 0 0;0 0 0 a b c];
K=(A*E/l)*T'*[1 -1;-1 1]* T;
end

