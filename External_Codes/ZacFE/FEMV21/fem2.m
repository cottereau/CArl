function [u,x,y,c]=fem2(nod2xy,el2nod, alpha,beta,s, bd)
%FEM2  Finite element solver in 2D.
%   FEM2 is a PDE solver using the Finite Element Method
%   for two dimensions.
%
%   U = FEM2(NOD2XY,EL2NOD,ALPHA,BETA,S,BD)
%            - for plotting with PLOT3(NOD2XY(:,1),NOD2XY(:,2),U,'.')
%
%   [U,X,Y,C] = FEM2(...)
%            - for plotting with FILL3(X,Y,U,C)
%
%   Here U is the solution, NOD2XY (#NODx2) is the coordinates of each node
%   (including boundary nodes) and EL2NOD (#ELx3) is the elements with
%   associated nodes
%   ALPHA, BETA and S are used in the PDE of the form:
%
%      -div(ALPHA*grad(U)) + BETA*U = S
%
%   Z is then satisfying the linear combination
%
%      U(x,y) = sum(Z_i*phi_i(x,y))
%
%   where phi_i is the i:th basis function (hat-function).
%   However the output U is Z.
%   U, ALPHA, BETA and S have lengths (#NOD).
%   BD is the boundary conditions entered as a cell array of
%   (#D_i)x2 or (#R_i)x3 -matrices as:
%
%      BD{i} = [N_1 DIR_1
%               ... ...
%               N_k DIR_k]
%   or
%
%      BD{i} = [N_1 ROB_1 GAMMA_1
%               ... ...   ...
%               N_k ROB_k GAMMA_k]
%
%   N_j (j=1..k=1..#D_i or 1..#R_i) is node number j
%   along boundary number i.
%   The integration is performed along the direction
%   of increasing row-indicies j.
%
%   The following shows how the constants are used:
%
%      Dirichlet:
%         U(i) = DIR_i
%      Robin:
%         n*(ALPHA*grad(U(i))) + GAMMA_i*U(i) = ROB_i
%
%   Neumann boundary condition is achieved by setting GAMMA_i
%   to zero.
%
%
%   See also GENMAT1, GENMAT2, FEM1, PLOTGRID2.

% Copyright (c) 2002-04-28, B. Rasmus Anthin.

error(nargchk(6,6,nargin))

if size(nod2xy,2)~=2
   error('NOD2XY must be a Nx2 matrix.')
elseif size(el2nod,2)~=3
   error('EL2NOD must be a Nx3 matrix.')
elseif prod(size(alpha))~=size(nod2xy,1)
   error('ALPHA must have as many elements as there are nodes.')
elseif prod(size(beta))~=size(nod2xy,1)
   error('BETA must have as many elements as there are nodes.')
elseif prod(size(s))~=size(nod2xy,1)
   error('S must have as many elements as there are nodes.')
end
alpha=alpha(:).';
beta=beta(:).';
s=s(:).';

[A2,B2,b]=genmat2(nod2xy,el2nod,alpha,beta,s);
A=A2+B2;

if length(bd)>=size(nod2xy,1), error('Too many boundary nodes. Solution already known?'),end

%Extract boundary conditions/nodes.
cond.d=[];
cond.r=cell(1);
j=1;
no.b=0;N.b=[];
no.r={};N.r={};
no.d=0;N.d=[];
if ~iscell(bd)
   tmp=cell(1);
   tmp{1}=bd;
   bd=tmp;
   clear tmp
end
for i=1:length(bd)
   no.b=no.b+size(bd{i},1);       %boundary nodes
   N.b=[N.b bd{i}(:,1)'];         %boundary nodes
   if size(bd{i},2)==2
      cond.d=[cond.d;bd{i}(:,2)];
      no.d=no.d+size(bd{i},1);    %Dirichlet nodes
      N.d=[N.d bd{i}(:,1)'];      %Dirichlet nodes
   elseif size(bd{i},2)==3
      cond.r{j}=bd{i}(:,2:3);
      no.r{j}=size(bd{i},1);      %Robin nodes
      N.r{j}=bd{i}(:,1)';         %Robin nodes
      j=j+1;
   else
      error('Boundary conditions must be a #Dx2 or #Rx3 matrix.')
   end
end

no.a=size(nod2xy,1);     %all nodes
N.a=1:no.a;              %all nodes

N.nd=setdiff(N.a,N.d);   %non Dirichlet nodes
no.nd=length(N.nd);      %non Dirichlet nodes

N.i=setdiff(N.a,N.b);    %interior nodes
no.i=no.a-no.b;          %interior nodes


C=A(N.nd,N.d)*sparse(cond.d);
if isempty(C),C=spalloc(N.nd,1,0);end
[D,E]=robin(cond.r,nod2xy,N,no);
A=A(N.nd,N.nd)+E(N.nd,N.nd);
b=b(N.nd)-C+D(N.nd);
clear C D E
z=A\b;
z=full(z);
u(N.nd)=z;
u(N.d)=cond.d;
if any(nargout==[3 4])
   I=el2nod';
   x=reshape(nod2xy(I,1),size(I));
   y=reshape(nod2xy(I,2),size(I));
   u=reshape(u(I),size(I));
   c=mean(u);
end



function [D,E]=robin(R,nod2xy,N,no)
D=spalloc(no.a,1,0);
E=spalloc(no.a,no.a,0);
if ~isempty(N.r)
   %Calculate pseudo-grid for each set of nodes
   %and calculate corresponding B_ij-matrixes
   for i=1:length(N.r)
      G{i}=cumsum([0 sqrt(sum(diff(nod2xy(N.r{i},:)',[],2).^2,1))]);
      [A,B]=genmat1(G{i},zeros(size(G{i})),ones(size(G{i})),zeros(size(G{i})));
      mat{i}=B;
   end
   %Calculate matrices D and E
   for i=1:length(N.r)
      D(N.r{i})=mat{i}*R{i}(:,1);
      E(N.r{i},N.r{i})=mat{i}.*R{i}(:,2*ones(1,length(R{i})))';
   end
end
