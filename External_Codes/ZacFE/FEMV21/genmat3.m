function [A,B,b]=genmat3(xyz,el,a,b,s)
%GENMAT3  Generates FEM3 matrices.
%   [A,B,b] = GENMAT3(NOD2XYZ,EL2NOD,ALPHA,BETA,S)
%   where NOD2XYZ contains the node points
%   and EL2NOD contains the tetrahedroidal elements (see FEM3).
%   ALPHA, BETA, S are parameters in the elliptic PDE:
%
%      -div(ALPHA*grad(U)) + BETA*U = S
%
%   A, B and b are matrices returned for solving
%   the PDE from the main routine FEM2 generally by using
%   the relation (A+B)*z=b.
%
%   See also FEM3.

% Copyright (c) 2002-04-14, B. Rasmus Anthin.

%A_ij = int(a*phi_i'*phi_j')dS
%B_ij = int(b*phi_i*phi_j)dS
%C_ij = int(phi_i*phi_j)dS
%A = A_ij
%B = B_ij
%b = C_ij*s'

aa=sparse(a(:).');
bb=sparse(b(:).');
ss=sparse(s(:).');
x=reshape(xyz(el,1),size(el))';
y=reshape(xyz(el,2),size(el))';
z=reshape(xyz(el,3),size(el))';
N=size(xyz,1);
A=spalloc(N,N,0);
B=spalloc(N,N,0);
b=spalloc(N,1,0);
for i=1:size(el,1)
   n=el(i,:);
   [Ae,Be]=cmpel(x(:,i),y(:,i),z(:,i));
   aaa=aa([1 1 1 1],n);
   bbb=bb([1 1 1 1],n);
   sss=ss(n);
   A(n,n)=A(n,n)+aaa.*Ae;
   B(n,n)=B(n,n)+bbb.*Be;
   b(n,1)=b(n,1)+Be*sss.';
end


function [Ae,Be]=cmpel(x,y,z)
V=elvol(x,y,z);
for i=1:4
   [a(i) b(i) c(i) d(i)]=cmpphi(x,y,z,i);
end
Ae=[b;c;d].'*[b;c;d]*V;
%Be=S/12*(eye(3)+1);
Be=intphi(x,y,z);
Ae=sparse(Ae); Be=sparse(Be);

function [a,b,c,d]=cmpphi(x,y,z,n)
V=elvol(x,y,z);
i=n+1; i=i-4*(i>4);
j=n+2; j=j-4*(j>4);
k=n+3; k=k-4*(k>4);
I=[i j k];
J=[j k i j];
x=x(:).';
y=y(:).';
z=z(:).';
dx=diff(x(J));
dy=diff(y(J));
dz=diff(z(J));
a=(-1)^i*det([x(I);y(I);z(I)])/3/V;
b=(-1)^i*y(I)*dz'/3/V;
c=(-1)^i*z(I)*dx'/3/V;
d=(-1)^i*x(I)*dy'/3/V;

% gives the same result as S/12*(eye(3)+1)
% this is the matrix for int(phi_i*phi_j)
function I=intphi(x,y,z)
V=elvol(x,y,z);
for i=1:4
   for j=1:4
      a=(i==1)+(j==1);
      b=(i==2)+(j==2);
      c=(i==3)+(j==3);
      d=(i==4)+(j==4);
      I(i,j)=3*V*prod(1:a)*prod(1:b)*prod(1:c)*prod(1:d)/prod(1:(a+b+c+d+3));
   end
end

function V=elvol(x,y,z)
r1=[x(1) y(1) z(1)];
r2=[x(2) y(2) z(2)];
r3=[x(3) y(3) z(3)];
r4=[x(4) y(4) z(4)];
V=(r1-r2)*(cross(r4-r2,r3-r2)).'/3;
