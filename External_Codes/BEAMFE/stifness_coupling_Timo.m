function [KC1,KC2] = stifness_coupling_Timo( elements,coordinates)


NE=size(elements,1); %number of elements
DIM=size(coordinates,2); %problem dimension
gaussLocations=[-1 1 -1];

dfunc = shapeder ([0 0;1 0;0 1]','P1');
dfunu =repmat([1;-1;1],1,3);%*invJ0
dfunt=dfunu;
funt=([1-gaussLocations;1+gaussLocations;1-gaussLocations]/2);


coordB = coordinates(coordinates(:,2)==0,1);
NEB = length(coordinates);
KC1 = zeros(2*length(coordinates),2*length(coordinates));
KC2 = zeros(length(coordinates),2*length(coordinates));
for ijk = 1:NE
    indice = elements(ijk,:);
    xelement = unique(coordinates(indice,1));
   % [~,indiceB] = ismember(xelement,coordB);%[find(xelement(1)-coordB==0,1) find(xelement(2)-coordB==0,1)];
   % indiceB(indiceB==0)=[];
    indiceM11 = [indice indice+NEB];
    indiceM12 = [2*indice-1 2*indice];
detJ0=(polyarea(coordinates(indice,1),coordinates(indice,2)));
mattemp = zeros(2*length(indice),2*length(indice));
invJ0 = sqrt(1/(max(xelement)-min(xelement)));
mattemp(1:3,[1 3 5]) =dfunu*squeeze(dfunc(1,:,:))';
mattemp(4:6,[2 4 6]) =dfunu*squeeze(dfunc(1,:,:))';
KC1(indiceM11,indiceM12) = KC1(indiceM11,indiceM12)+mattemp;


indiceM21 = [indice];
indiceM22 = [2*indice-1];
X2 = repmat((coordinates(indice,2))',3,1);
mattemp1 =-X2.*dfunt*squeeze(dfunc(1,:,:))';
mattemp2 =-funt*squeeze(dfunc(2,:,:))'*detJ0;
KC2(indice,indiceM22) = mattemp1;
indiceM22 = [2*indice];
KC2(indice,indiceM22) = KC2(indice,indiceM22)+mattemp2;
end