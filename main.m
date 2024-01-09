% main script 
coordinates=load('coordi.txt');
connectivity=load('connectivity.txt');

% known forces
kf=load('kf.txt');

% known displacements
kd=load('ku.txt');

%number of nodes
a=size(coordinates); 

% dof is thrice of number of nodes
m=3*a(1,1);

%size of global stiffness matrix denoted by n
n=m;

K_g=zeros(m,m);
a=size(connectivity);
m=a(1,1);
for i=1:m
  node_1=connectivity(i,1);
  node_2=connectivity(i,2);
  x1=coordinates(node_1,1);
  y1=coordinates(node_1,2);
  z1=coordinates(node_1,3);
  x2=coordinates(node_2,1);
  y2=coordinates(node_2,2);
  z2=coordinates(node_2,3);
  A=200*1e-6;
  E=70*1e6;
  k_e = ElementStiffness(x1,y1,z1,x2,y2,z2,A,E);
  %key=[1 2 3 4 5 6];
 % value=[3*node_1-2 3*node_1-1 3*node_1 3*node_2-2 3*node_2-1 3*node_2];
  %M=containers.Map(key,value);
%for j=1:6
  %    for k=1:6
  %  K_g(M(j),M(k))=K_g(M(j),M(k))+k_e(j,k);
  %    end
  %end 
  
  ig=[3*node_1-2:3*node_1,3*node_2-2:3*node_2];
  jg=ig;
  K_g(ig,jg)=K_g(ig,jg)+k_e(1:6,1:6);
end
% construct matrices for solution 
% defining two matrix c,d for lower division of K_g
a=size(kf);
b=size(kd);
p=a(1,1);
q=b(1,1);
C=zeros(p,n-p);
D=zeros(p,p);

% construction of C
for i=1:p
    for j=1:n-p
        C(i,j)=K_g(kf(i,1),kd(j,1));
    end
end

% construction of D
for i=1:p
    for j=1:p
        D(i,j)=K_g(kf(i,1),kf(j,1));
    end
end

f=zeros(p,1);
for i=1:p
    f(i,1)=kf(i,2);
end
ku=zeros(q,1);
for i=1:q
    ku(i,1)=kd(i,2);
end
f=f-C*(ku);
D=inv(D);

% unknown displacements
d=D*f;

% displacement vector
disp=zeros(n,1);
i=1;
j=1;
for k=1:n
    if i<=q && kd(i,1)==k
        disp(k,1)=kd(i,2);
        i=i+1;
    else 
        disp(k,1)=d(j,1);
        j=j+1;
    end 
end


% force matrix
force=K_g*disp;


% K_g
fprintf('The outputs for displacement and force are printed according to their serial degree of freedoms\n');


% displacement at node is stored in the variable disp
disp


% forces at any dof is stored in the variable force
force

%Plot of undeformed structure 
s=size(connectivity);
e=s(1,1);
figure
f=10;%factor to enlarge the graph
for i=1:e
    node_1=connectivity(i,1);
    node_2=connectivity(i,2);
    x=[coordinates(node_1,1),coordinates(node_2,1)];
    y=[coordinates(node_1,2),coordinates(node_2,2)];
    z=[coordinates(node_1,3),coordinates(node_2,3)];
    plot3(x,y,z,'-ob','LineWidth',1);
    hold on 
end
plot3(0,-5,0);
%Plot of deformed structure 
for i=1:e
    node_1=connectivity(i,1);
    node_2=connectivity(i,2);
    x=[coordinates(node_1,1)+disp((node_1-1)*3+1,1)*f,coordinates(node_2,1)+disp((node_2-1)*3+1,1)*f];
    y=[coordinates(node_1,2)+disp((node_1-1)*3+2,1)*f,coordinates(node_2,2)+disp((node_2-1)*3+2,1)*f];
    z=[coordinates(node_1,3)+disp((node_1-1)*3+3,1)*f,coordinates(node_2,3)+disp((node_2-1)*3+3,1)*f];
    plot3(x,y,z,'-or','LineWidth',1); 
    hold on
end
title('UNDEFORMED SHAPE DEFORMED SHAPE');
hold off
