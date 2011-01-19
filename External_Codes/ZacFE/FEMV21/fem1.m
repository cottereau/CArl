function u=fem1(xn, alpha,beta,s, bd)
%FEM1  Finite element solver in 1D.
%   FEM1 is a PDE solver using the Finite Element Method
%   for one dimension.
%
%   U = FEM1(XN,ALPHA,BETA,S,BD)
%
%   where U is the solution and XN is the gridpoints (including endpoints)
%   for the mesh.
%   ALPHA, BETA and S are used in the PDE of the form:
%
%      -(ALPHA*U')' + BETA*U = S
%
%   Z is then satisfying the linear combination
%
%      U(x) = sum(Z_i*phi_i(x))
%
%   where phi_i is the i:th basis function (hat-function).
%   However the output U is Z.
%   BD is the boundary conditions entered as a 2x3-matrix as:
%
%      BD = [DIRa ROBa GAMMAa
%            DIRb ROBb GAMMAb]
%
%   Any type of condition not used have to be replaced by a NaN.
%   That is if using Dirichlet for endpoint "a" and Neumann for
%   endpoint "b" you can write:
%
%      BD = [1.5 nan nan
%            nan -10   0]
%
%   The following shows how the constants are used:
%
%      Dirichlet:
%         U(a) = DIRa
%         U(b) = DIRb
%      Robin:
%         U'(a) + GAMMAa*U(a) = ROBa
%         U'(b) + GAMMAb*U(b) = ROBb
%
%   Neumann boundary condition is achieved by setting GAMMAa or
%   GAMMAb to zero.
%
%
%   See also GENMAT1, REFINE1.

% Copyright (c) 2002-03-13, B. Rasmus Anthin.

error(nargchk(5,5,nargin))

xn=xn(:).';
alpha=alpha(:).';
beta=beta(:).';
s=s(:).';
N=length(xn);

%old scheme:
%-----------
%for i=1:N
%   for j=1:N
%      A(i,j)=trapz(x,alpha.*dphi(i,:).*dphi(j,:)+beta.*phi(i,:).*phi(j,:));
%   end
%end
%for i=1:N
%   b(i,1)=trapz(xn,phi(i,:).*s);
%end

[A,B,b]=genmat1(xn,alpha,beta,s);
A=A+B;

if size(bd)~=[2 3], error('Boundary conditions must be a 2x3 matrix.'),end

isbd=~isnan(bd);
a1=alpha(1); a2=alpha(end);
b1=beta(1); b2=beta(end);
[i j]=ind2sub(size(bd),find(isbd));
bd2=bd(sub2ind(size(bd),i(j~=3),j(j~=3)));
gam=bd(:,3);

if sum(isbd(1,1:2))~=1 | sum(isbd(2,1:2))~=1
   error('Must enter exactly one boundary condition on each endpoint.')
elseif isbd(1,2) & ~isbd(1,3)
   error('Must define GAMMAa for Robin/Neumann condition.')
elseif isbd(2,2) & ~isbd(2,3)
   error('Must define GAMMAb for Robin/Neumann condition.')
end

M=N-length(find(isbd(:,1)));
ii=1+isbd(1,1) : N-isbd(2,1);
if isbd(1,1), kk=1; else, kk=[]; end
if isbd(2,1), kk=[kk N]; end
mm=find(isbd(:,1));
C=A(ii,kk)*sparse(bd(mm,1));
D=spalloc(M,1,2);
if isbd(1,2), D(1)=-bd2(1)*a1; end
if isbd(2,2), D(M)=bd2(2)*a2; end
E=spalloc(M,M,2);
if isbd(1,2), E(1)=-a1*gam(1);end
if isbd(2,2), E(M,M)=a2*gam(2);end
A=A(ii,ii)+E;
b=b(ii)-C+D;
clear C D E ii
z=A\b;
z=full(z);
if isbd(1,1), z=[bd2(1);z]; end
if isbd(2,1), z=[z;bd2(2)]; end
u=z;