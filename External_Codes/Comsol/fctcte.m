function [upc,in] = fctcte(Xf, Tf, uf, pts ,pcn)
% function [upc,in] = fctintp(Xf, Tf, uf, pts ,pcn)
% function compute the linear interpolation for upc from uf and from the closer nodes
% nodes (in Xf) of the element Tf where the point pc is.
% system
% ------------------------------------------------------------------------
% Input variables:
%           Xf   		 : nodes of the fine mesh
%           Tf   		 : connectivity of Xf
%           uf    		 : solution on the fine mesh
%           pts           : vector of any points (np*2)
%           pcn          : value of upc if pc is not in an element of Tf
% Output variables:
%           upc    		 : interpoled solution
%           in           : in = 0 if pc is not in an element of Tf
% Files loaded
%           none
% Files generated 
%           none
%--------------------------------------------------------------------------

% Cedric Zaccardi
% cedric.zaccardi@ecp.fr
% Ecole Centrale Paris
% LMSSMAT Laboratoire de Mecanique de Sols, Structures et Materiaux
% Grande Voie des Vignes
% F-92295 Chatenay-Malabry Cedex
% France

[Np,d] = size(pts) ;
if d~=2
    error('Wrong dimension for the use of interp2d')
end

in = zeros(Np,1) ;
upc = zeros(Np,1) ;

for i2 = 1:Np
    pcx = pts(i2,1) ;
    pcy = pts(i2,2) ;
    indi = zeros(size(Tf,1),1) ;

    for i1 = 1:size(Tf,1)
        a = Tf(i1,1) ;
        b = Tf(i1,2) ;
        c = Tf(i1,3) ;
        pXkx = [Xf(a,1) ; Xf(b,1) ; Xf(c,1) ; Xf(a,1)] ;
        pXky = [Xf(a,2) ; Xf(b,2) ; Xf(c,2) ; Xf(a,2)] ;
        
         % remarque sur le 0.0003 : c'est pour assurer le bon fonctionnement de
         % inpolygon, mais le parametre est a regler
         in1 = inpolygon(pcx+0.003 , pcy , pXkx , pXky) ;
         in2 = inpolygon(pcx-0.003 , pcy , pXkx , pXky) ;
         in3 = inpolygon(pcx , pcy+0.003 , pXkx , pXky) ;
         in4 = inpolygon(pcx , pcy-0.003 , pXkx , pXky) ;
         indi(i1,1) = in1 | in2 | in3 | in4 ;
    end

    clear a b c ;

    ind = find(indi) ;

    if size(ind,1)~=0
        ind = ind(1,1) ;
        in(i2,1) = 1 ;

        a = Tf(ind,1) ;
        b = Tf(ind,2) ;
        c = Tf(ind,3) ;

        upc(i2,1) = (uf(a,1) + uf(b,1) + uf(c,1)) / 3  ;
    else
        in(i2,1) = 0 ;
        upc(i2,1) = pcn ;
    end

end