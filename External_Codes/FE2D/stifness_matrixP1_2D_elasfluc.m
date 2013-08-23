function [A,areas]=stifness_matrixP1_2D_elasfluc(elements,coordinates,lambda,mu,alpha)

%for laplace
NE=size(elements,1); %number of elements
DIM=size(coordinates,2); %problem dimension
 
%laplace 
NLB=3; %number of local basic functions, it must be known!
coord=zeros(DIM,NLB,NE);
for d=1:DIM
     for i=1:NLB
         coord(d,i,:)=coordinates(elements(:,i),d);
     end
end
%%%%%%%%%%%%%%%%%%%
% IP = [1/3;1/3];
% WG = 1/2;
%%%%%%%%%%%%%%%%%%%
IP = [1/6 2/3 1/6; 2/3 1/6 1/6];
WG = 1/6;
%%%%%%%%%%%%%%%%%%%
NG = size(IP,2); % number of gauss points
[dphi,jac] = phider(coord,IP,'P1'); %integration rule, it must be known!
jac = abs(squeeze(jac))'.* alpha;
areas = reshape( jac, 1, NG*NE )*WG;

%elasticity matrix is derived from laplace matrix
%R = zeros(3,6,NE);
dphi = reshape( permute(dphi,[1 2 4 3]), DIM, NLB, NG*NE );
R = zeros( NLB, NLB*DIM, NG*NE );
R([1,3],[1,3,5],:) = dphi;  
R([3,2],[2,4,6],:) = dphi;
clear dphi

if numel(lambda)~=1
    lambda = (lambda(:,:,1) * [1/6 2/3 1/6; 2/3 1/6 1/6;1/6 1/6 2/3])';
else
    lambda = lambda*ones(NG,NE);
end
if numel(mu)~=1
    mu = (mu * [1/6 2/3 1/6; 2/3 1/6 1/6;1/6 1/6 2/3])';
else
    mu = mu*ones(NG,NE);
end
C = mu(:)*[2 0 0 0 2 0 0 0 1] + lambda(:)*[1 1 0 1 1 0 0 0 0];
C = reshape( C', 3, 3, NG*NE );

Elements=2*elements(:,[1 1 2 2 3 3])-kron(ones(NE,1),[1,0,1,0,1,0]);
Y=reshape(repmat(Elements,1,NLB*DIM)',NLB*DIM,NLB*DIM,NE);
Y = repmat( Y, [1 1 NG] );
X=permute(Y,[2 1 3]);
Z=astam(areas',amtam(R,amtam(C,R)));
% Z=astam(areas',amtam(R,smamt(C,permute(R,[2 1 3]))));

A=sparse(X(:),Y(:),Z(:));

end


