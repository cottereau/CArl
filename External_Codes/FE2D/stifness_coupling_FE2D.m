function [KC1,KC2] = stifness_coupling_FE2D( elements,coordinates)


NE=size(elements,1); %number of elements
DIM=size(coordinates,2); %problem dimension
gaussLocations=[1/3];

dfunc = shapeder ([1/3;1/3],'P1');
dfunu =[1;-1;-1];%*invJ0
dfunt=dfunu;
funt=([gaussLocations;1-gaussLocations;1-gaussLocations]);


%coordB = coordinates(coordinates(:,2)==0,1);
NEB = length(coordinates);
KC1 = zeros(2*length(coordinates),2*length(coordinates));
KC2 = zeros(2*length(coordinates),length(coordinates));
for ijk = 1:NE
    indice = elements(ijk,:);
    xelement = sqrt((coordinates(indice,1)));
        dist = sqrt(abs(diff(xelement)));
    dist(end+1)=sqrt(abs(xelement(end)-xelement(1)));
   % [~,indiceB] = ismember(xelement,coordB);%[find(xelement(1)-coordB==0,1) find(xelement(2)-coordB==0,1)];
   % indiceB(indiceB==0)=[];
    indiceM11 = [indice indice+NEB];
    indiceM12 = [2*indice-1 2*indice];
detJ0=(polyarea(coordinates(indice,1),coordinates(indice,2)));
mattemp = zeros(2*length(indice),2*length(indice));
%invJ0 = sqrt(1/(max(xelement)-min(xelement)));
mattemp([1 3 5],1:3) =kron(squeeze(dfunc(1,:,:))',dfunu');
mattemp([2 4 6],4:6) =kron(squeeze(dfunc(1,:,:))',dfunu');
KC1(indiceM12,indiceM11) = KC1(indiceM12,indiceM11)+mattemp;


%indiceM21 = [indice];
indiceM22 = [2*indice-1];
X2 = coordinates(indice,2);
mattemp1 =-kron(squeeze(dfunc(1,:,:))',(X2.*dfunt)');
mattemp2 =-kron(squeeze(dfunc(2,:,:))',(funt.*dist)');
KC2(indiceM22,indice) = KC2(indiceM22,indice) + mattemp1;
indiceM22 = [2*indice];
KC2(indiceM22,indice) = KC2(indiceM22,indice)+mattemp2;
end