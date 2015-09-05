function [A,areasave]=stifness_matrixP1_2D_cubic(elements,coordinates,c1,c4,c12,alpha,theta)

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
areasave = abs(squeeze(jac))/factorial(DIM);
areasave = areasave(1,:)';
if size(jac,3)==1
    jac = abs(squeeze(jac)).* alpha;
else
    jac = abs(squeeze(jac))'.* alpha;
end
areas = reshape( jac, 1, NG*NE )*WG;

%elasticity matrix is derived from laplace matrix
%R = zeros(3,6,NE);
dphi = reshape( permute(dphi,[1 2 4 3]), DIM, NLB, NG*NE );
R = zeros( NLB, NLB*DIM, NG*NE );
R([1,3],[1,3,5],:) = dphi;  
R([3,2],[2,4,6],:) = dphi;
clear dphi

if numel(c1)~=1
    c1 = (c1 * [1/6 2/3 1/6; 2/3 1/6 1/6;1/6 1/6 2/3])';
else
    c1 = c1*ones(NG,NE);
end
c4 = c4*ones(NG,NE);
c12 = c12*ones(NG,NE);
if nargin>6
    ct = ones(NG,1)*cos(theta(:))';
    st = ones(NG,1)*sin(theta(:))';
else
    ct = ones(NG,NE);
    st = zeros(NG,NE);
end

% C = [c1(:) c12(:) zeros(NG*NE,1) c12(:) c1(:) zeros(NG*NE,3) c4(:)];
c4ps4 = ct(:).^4+st(:).^4;
c2s2 = ct(:).^2.*st(:).^2;
c2ms2cs = (ct(:).^2-st(:).^2).*ct(:).*st(:);
C = [c1(:).*c4ps4+2*c12(:).*c2s2+4*c4(:).*c2s2 ...
     2*c1(:).*c2s2+c12(:).*c4ps4-4*c4(:).*c2s2 ...
     -c1(:).*c2ms2cs+c12(:).*c2ms2cs+2*c4(:).*c2ms2cs ...
     2*c1(:).*c2s2+c12(:).*c4ps4-4*c4(:).*c2s2 ...
     c1(:).*c4ps4+2*c12(:).*c2s2+4*c4(:).*c2s2 ...
     c1(:).*c2ms2cs-c12(:).*c2ms2cs-2*c4(:).*c2ms2cs ...
     -c1(:).*c2ms2cs+c12(:).*c2ms2cs+2*c4(:).*c2ms2cs ...
     c1(:).*c2ms2cs-c12(:).*c2ms2cs-2*c4(:).*c2ms2cs ...
     2*c1(:).*c2s2-2*c12(:).*c2s2+c4(:).*(ct(:).^2-st(:).^2).^2];
C = reshape( C', 3, 3, NG*NE );

Elements=2*elements(:,[1 1 2 2 3 3])-kron(ones(NE,1),[1,0,1,0,1,0]);
Y=reshape(repmat(Elements,1,NLB*DIM)',NLB*DIM,NLB*DIM,NE);
Y = repmat( Y, [1 1 NG] );
X=permute(Y,[2 1 3]);
Z=astam(areas',amtam(R,amtam(C,R)));
% Z=astam(areas',amtam(R,smamt(C,permute(R,[2 1 3]))));

A=sparse(X(:),Y(:),Z(:));

end


