function [ M1,M2 ] = mass_coupling_FE2D(elements,coordinates)

NE=size(elements,1); %number of elements
DIM=size(coordinates,2); %problem dimension
func = shapefun([1/3;1/3],'P1');
gaussLocations=[1/3];
funu = ([gaussLocations;1-gaussLocations;1-gaussLocations]);   
funt=funu;

%coordB = (coordinates(find(coordinates(:,2)==0),1));
coordB = coordinates(coordinates(:,2)==0,1);
NEB = length(coordinates);
M1 = zeros(2*length(coordinates),2*length(coordinates));
M2 = zeros(2*length(coordinates),length(coordinates));
for ijk = 1:NE
    indice = elements(ijk,:);
    xelement = (coordinates(indice,1));
    dist = sqrt(abs(diff(xelement)));
    dist(end+1)=sqrt(abs(xelement(end)-xelement(1)));
 %   [~,indiceB] = ismember(xelement,coordB);%[find(xelement(1)-coordB==0,1) find(xelement(2)-coordB==0,1)];
  %  indiceB(indiceB==0)=[];
    indiceM11 = [indice indice+NEB];
    indiceM12 = [2*indice-1 2*indice];
detJ0=(polyarea(coordinates(indice,1),coordinates(indice,2)));
mattemp = zeros(2*length(indice),2*length(indice));
mattemp([1 3 5],1:3) =kron(func,(funu.*dist)')*sqrt(detJ0);
mattemp([2 4 6],4:6) =kron(func,(funu.*dist)')*sqrt(detJ0);
M1(indiceM12,indiceM11) = M1(indiceM12,indiceM11)+mattemp;


%indiceM21 = [indice];
indiceM22 = [2*indice-1];
X2 = coordinates(indice,2);
mattemp =-kron(func,(X2.*(funt.*dist))')*sqrt(detJ0);
M2(indiceM22,indice) = M2(indiceM22,indice)+mattemp;




end