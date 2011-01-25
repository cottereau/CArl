function [Ke,fe] = elementStiffnessMatrixHomeFE( Xe, Epg, numberOfNodes, ...
                           pospg, pespg, N, Nxi, Neta, load ) 
% [Ke,fe] = elementStiffnessMatrixHomeFE( Xe, numberOfNodes, pospg, pespg, ...
%   N, Nxi, Neta ) 
%
% creates an elemental matrix
%
% INPUT
%   Xe             nodal coords
%   numberOfNodes  number of element nodes
%   pospg,pespg   position and weigth of gauss points 
%   N,Nxi,Neta    shape functions and its derivetives
%
% OUTPUT
%   K
%   f
%

d = size(Xe,2);
numberOfGaussPoints = size( pospg, 1 );

Ke = zeros( numberOfNodes, numberOfNodes );
fe = zeros( numberOfNodes, size(load,2) );

if d==1
    for igaus = 1:numberOfGaussPoints
        dN = Nxi(igaus,:);
        jacob = dN*Xe ;
        dvolu = pespg(igaus) * det( jacob );
        res = jacob \ dN;
        Nx = res(1,:);
        Ke = Ke + ( Nx'*Nx ) * dvolu * Epg(igaus);
        fe = fe + N(igaus,:)' * load(igaus,:) * dvolu;
    end

elseif d==2
    for igaus = 1:numberOfGaussPoints
        dN = [ Nxi(igaus,:) ; Neta(igaus,:) ];
        jacob = dN*Xe ;
        dvolu = pespg(igaus) * det( jacob );
        res = jacob \ dN;
        Nx = res(1,:);
        Ny = res(2,:);
        Ke = Ke + ( Nx'*Nx + Ny'*Ny ) * dvolu * Epg(igaus);
        fe = fe + N(igaus,:)' * load(igaus,:) * dvolu;
    end
end

