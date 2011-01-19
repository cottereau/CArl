function [upc,in] = fctintp(Xf, Tf, uf, pts ,pcn)
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
    pc = pts(i2,:) ;
    
    indi = zeros(size(Tf,1),1) ;

    for i1 = 1:2 ;%size(Tf,1)
        a = Tf(i1,1) ;
        b = Tf(i1,2) ;
        c = Tf(i1,3) ;
        pXk = [Xf(a,:) ; Xf(b,:) ; Xf(c,:)] ;

         % remarque sur le 0.0003 : c'est pour assurer le bon fonctionnement de
         % inpolygon, mais le parametre est a regler

         in1 = inpolygon(pc(:,1)+0.003,pc(:,2),pXk(:,1),pXk(:,2)) ;
         in2 = inpolygon(pc(:,1)-0.003,pc(:,2),pXk(:,1),pXk(:,2)) ;
         in3 = inpolygon(pc(:,1),pc(:,2)+0.003,pXk(:,1),pXk(:,2)) ;
         in4 = inpolygon(pc(:,1),pc(:,2)-0.003,pXk(:,1),pXk(:,2)) ;
         indi(i1,1) = in1 | in2 | in3 | in4 ;
    end

    clear a b c ;

    ind = find(indi) ;
    ind = ind(1,1) ;

    if size(ind,1)==1

        in(i2,1) = 1 ;

        a = Tf(ind,1) ;
        b = Tf(ind,2) ;
        c = Tf(ind,3) ;

        Vab = Xf(b,:) - Xf(a,:) ;
        Vac = Xf(c,:) - Xf(a,:) ;
        Vax = pc - Xf(a,:) ;

        detr = (Vab(1,1)*Vac(1,2) - Vab(1,2)*Vac(1,1)) ;
        Alpha = (1/detr)*[Vac(1,2),-Vac(1,1);-Vab(1,2),Vab(1,1)]...
             *Vax(1,:)' ;
        alphaab = Alpha(1,:) ;
        alphaac = Alpha(2,:) ;

        upc(i2,1) = uf(a,1) + alphaab * ( uf(b,1) - uf(a,1)) + ...
                 alphaac * (uf(c,1) - uf(a,1)) ;
    else
        in(i2,1) = 0 ;
        upc(i2,1) = pcn ;
    end

end