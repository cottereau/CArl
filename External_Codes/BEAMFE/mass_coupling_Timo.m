function [ M1,M2 ] = mass_coupling_Timo(elements,coordinates)

NE=size(elements,1); %number of elements
DIM=size(coordinates,2); %problem dimension
func = shapefun([0 0;1 0;0 1]','P1');
gaussLocations=[-1 1 -1];
funu = ([1-gaussLocations;1+gaussLocations;1-gaussLocations]/2);   
funt=funu;

%coordB = (coordinates(find(coordinates(:,2)==0),1));
coordB = coordinates(coordinates(:,2)==0,1);
NEB = length(coordinates);
M1 = zeros(2*length(coordinates),2*length(coordinates));
M2 = zeros(length(coordinates),2*length(coordinates));
for ijk = 1:NE
    indice = elements(ijk,:);
    xelement = unique(coordinates(indice,1));
 %   [~,indiceB] = ismember(xelement,coordB);%[find(xelement(1)-coordB==0,1) find(xelement(2)-coordB==0,1)];
  %  indiceB(indiceB==0)=[];
    indiceM11 = [indice indice+NEB];
    indiceM12 = [2*indice-1 2*indice];
detJ0=(polyarea(coordinates(indice,1),coordinates(indice,2)));
mattemp = zeros(2*length(indice),2*length(indice));
mattemp(1:3,[1 3 5]) =funu*func'*detJ0;
mattemp(4:6,[2 4 6]) =funu*func'*detJ0;
M1(indiceM11,indiceM12) = M1(indiceM11,indiceM12)+mattemp;


indiceM21 = [indice];
indiceM22 = [2*indice-1];
X2 = repmat((coordinates(indice,2))',3,1);
mattemp =-X2.*funt*(func)'*detJ0;
M2(indice,indiceM22) = M2(indice,indiceM22)+mattemp;




end