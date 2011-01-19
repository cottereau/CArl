function [A,B,b]=genmat2(xy,el,a,b,s)
%GENMAT2  Generates FEM2 matrices.
%   [A,B,b] = GENMAT2(NOD2XY,EL2NOD,ALPHA,BETA,S)
%   where NOD2XY contains the node points
%   and EL2NOD contains the triangular elements (see FEM2).
%   ALPHA, BETA, S are parameters in the elliptic PDE:
%
%      -div(ALPHA*grad(U)) + BETA*U = S
%
%   A, B and b are matrices returned for solving
%   the PDE from the main routine FEM2 generally by using
%   the relation (A+B)*z=b.
%
%   See also FEM2.

% Copyright (c) 2002-04-14, B. Rasmus Anthin.
% modified by C. Zaccardi 05/2010

%A_ij = int(a*phi_i'*phi_j')dS
%B_ij = int(b*phi_i*phi_j)dS
%C_ij = int(phi_i*phi_j)dS
%A = A_ij
%B = B_ij
%b = C_ij*s'

aa=sparse(a(:).');
bb=sparse(b(:).');
ss=sparse(s(:).');
x=reshape(xy(el,1),size(el))';
y=reshape(xy(el,2),size(el))';
N=size(xy,1);
A=spalloc(N,N,0);
B=spalloc(N,N,0);
b=spalloc(N,1,0);
[Nel nne ] = size(el);
for i=1:Nel
   n=el(i,:);
   [Ae,Be]=cmpel(x(:,i),y(:,i));
%    aaa=aa([1 1 1],n);
   aaa=aa(1,n);
%    bbb=bb([1 1 1],n);
   bbb=bb(1,n);
   sss=ss(n);
   % taking the mean a for the element i
   aaa = mean(full(aaa)) * ones(nne) ;
   bbb = mean(full(bbb)) * ones(nne) ;
   % not checked for s...
   A(n,n)=A(n,n)+aaa.*Ae;
   B(n,n)=B(n,n)+bbb.*Be;
   b(n,1)=b(n,1)+Be*sss.';
end


function [Ae,Be]=cmpel(x,y)
S=polyarea(x,y);
for i=1:3
   [a(i) b(i) c(i)]=cmpphi(x,y,i);
end
% keyboard
z=zeros(1,3);
Ae=[b;c;z].'*[b;c;z]*S;
Be=S/12*(eye(3)+1);
%Be=intphi(x,y);
Ae=sparse(Ae); Be=sparse(Be);

function [a,b,c]=cmpphi(x,y,n)
S=polyarea(x,y);
i=n+1; i=i-3*(i>3);
j=n+2; j=j-3*(j>3);
a=(x(i)*y(j)-x(j)*y(i))/2/S;
b=diff(y([i j]))/2/S;
c=diff(x([j i]))/2/S;

% gives the same result as S/12*(eye(3)+1)
% this is the matrix for int(phi_i*phi_j)
function I=intphi(x,y)
S=polyarea(x,y);
for i=1:3
   for j=1:3
      a=(i==1)+(j==1);
      b=(i==2)+(j==2);
      c=(i==3)+(j==3);
      I(i,j)=2*S*prod(1:a)*prod(1:b)*prod(1:c)/prod(1:(a+b+c+2));
   end
end
