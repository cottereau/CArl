function [A,B,b]=genmat1(xn,a,b,s)
%GENMAT1  Generates FEM1 matrices.
%   [A,B,b] = GENMAT1(XN,ALPHA,BETA,S) where XN is the,
%   grid-points and ALPHA, BETA, S are parameters in the
%   elliptic PDE:
%
%      -(ALPHA*U')' + BETA*U = S
%
%   A, B and b are matrices returned for solving
%   the PDE from the main routine FEM1 generally by using
%   the relation (A+B)*z=b.
%
%   See also FEM1, REFINE1.

% Copyright (c) 2002-03-13, B. Rasmus Anthin.
% modified by C. Zaccardi 05/2010

%A_ij = int(a*phi_i'*phi_j')dx
%B_ij = int(b*phi_i*phi_j)dx
%C_ij = int(phi_i*phi_j)dx
%A = A_ij
%B = B_ij
%b = C_ij*S'

xn=xn(:)';
N=length(xn);
dxn=diff(xn);
s=sparse(s);

%a_{i,i} = a*[1/(x_n-x_{n-1}) + 1/(x_{n+1}-x_n]
aii=[a(1)/dxn(1), ...
   a(2:end-1)./dxn(1:end-1)+a(3:end)./dxn(2:end), ...
   a(end)/dxn(end)];
%b_{i,i} = b*[x_{n+1}-x_{n-1}]/3
bii=[b(1)*dxn(1), ...
   b(2:end-1).*(xn(3:end)-xn(1:end-2)), ...
   b(end)*dxn(end)]/3;
%c_{i,i} = [x_{n+1}-x_{n-1}]/3
cii=[dxn(1), ...
   (xn(3:end)-xn(1:end-2)), ...
   dxn(end)]/3;
%a_{i,i+1} = -a/(x_{n+1}-x_n)
aii1=-a(2:end)./dxn;
%a_{i+1,i} = -a/(x_{n+1}-x_n)
% ---------------------------------modification
% ai1i=-a(1:end-1)./dxn;
ai1i=-a(2:end)./dxn;
% ---------------------------------end modification
%b_{i,i+1} = b*(x_{n+1}-x_n)/6
bii1=b(2:end).*dxn/6;
%b_{i+1,i} = b*(x_{n+1}-x_n)/6
bi1i=b(1:end-1).*dxn/6;
%c_{i,i+1} = (x_{n+1}-x_n)/6
cii1=dxn/6;
%c_{i+1,i} = (x_{n+1}-x_n)/6
ci1i=dxn/6;

%A=diag(aii) + diag(aii1,1) + diag(ai1i,-1);
A=spdiags([[ai1i 0].' aii.' [0 aii1].'],-1:1,N,N);
%B=diag(bii) + diag(bii1,1) + diag(bi1i,-1);
B=spdiags([[bi1i 0].' bii.' [0 bii1].'],-1:1,N,N);
%C=diag(cii) + diag(cii1,1) + diag(ci1i,-1);
C=spdiags([[ci1i 0].' cii.' [0 cii1].'],-1:1,N,N);
% keyboard
b=C*s.';