%This is an electrostatic example. It is a resistive element
%with a potential at the bottom equal to V=0 and at the top
%the electric potential is V=10. There is no net electric field
%towards the other sides of the square plate.

N=10;
%generate grid
xn=linspace(0,1,N);
[X,Y]=meshgrid(xn);
nod2xy=[X(:) Y(:)];
%generate triangulation/mesh
el2nod=delaunay(X(:),Y(:));
%find boundary nodes
geom.a=find(nod2xy(:,2)==0);      %lower boundary
geom.b=find(nod2xy(:,1)==1);      %right boundary
geom.c=find(nod2xy(:,2)==1);      %upper boundary
geom.d=find(nod2xy(:,1)==0);      %left boundary
geom.b=setdiff(geom.b,[geom.a;geom.c]);
geom.d=setdiff(geom.d,[geom.a;geom.c]);

%drawing mesh and gridpoints
plotgrid2(nod2xy,el2nod)
%plotting boundary nodes
figure
hold on
plot(nod2xy(geom.a,1),nod2xy(geom.a,2),'.-b')
plot(nod2xy(geom.b,1),nod2xy(geom.b,2),'.-r')
plot(nod2xy(geom.c,1),nod2xy(geom.c,2),'.-b')
plot(nod2xy(geom.d,1),nod2xy(geom.d,2),'.-r')
legend('Dirichlet','Neumann',0)

%generate boundary conditions
bd={};
bd{1}=[geom.a zeros(N,1)];               %lower boundary (Dirichlet)
bd{2}=[geom.b zeros(N-2,2)];             %right boundary (Neumann)
bd{3}=[geom.c 10*ones(N,1)];             %upper boundary (Dirichlet)
bd{4}=[geom.d zeros(N-2,2)];             %left boundary (Neumann)

%set PDE parameter values
alpha=ones(size(nod2xy,1),1);
beta=zeros(size(nod2xy,1),1);
s=ones(size(nod2xy,1),1);

u=fem2(nod2xy,el2nod,alpha,beta,s,bd);
%plot solution as points
figure
plot3(nod2xy(:,1),nod2xy(:,2),u,'.')
rotate3d
%plot solution as a surface
[u x y c]=fem2(nod2xy,el2nod,alpha,beta,s,bd);
figure
fill3(x,y,u,c)
shading flat
colormap summer
rotate3d
