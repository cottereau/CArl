function M=mass_matrixP1_2D_elasticity(elements,coordinates)

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
IP=[1/3;1/3];
[~,jac] = phider(coord,IP,'P1'); %integration rule, it must be known! 
areas=abs(squeeze(jac))/factorial(DIM);

elements_elasticity=2*elements(:,[1 1 2 2 3 3])-kron(ones(NE,1),[1,0,1,0,1,0]);
NLB=6;
Y_3D=reshape(repmat(elements_elasticity,1,NLB)',NLB,NLB,NE);
X_3D=permute(Y_3D,[2 1 3]);

Zlocal=zeros(6);
Zlocal([1 3 5],[1 3 5])=(ones(3)+eye(3))/12;
Zlocal([2 4 6],[2 4 6])=(ones(3)+eye(3))/12;
Z_3D=astam(areas',repmat(Zlocal,[1 1 NE]));

M=sparse(X_3D(:),Y_3D(:),Z_3D(:));