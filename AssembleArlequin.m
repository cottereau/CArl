function [ K, F, opt ] = AssembleArlequin( model, coupling )
% ASSEMBLEARLEQUIN to assemble the Arlequin system before resolution
%
%  syntax: [ K, F, ind, opt ] = AssembleArlequin( Ki, Fi, C1, C2, ...
%                                                    model, coupling )
%
% the format is sparse so all matrices of Stiffness, Force, coupling are
% given in a (x,y,val) triplet (see SPARSE)
% All matrices are also cells, for the possibility of having several models
% and several couplings
%
% in opt, the beginning and ending indices are kept for each part of the
% matrix: K indicates the main matrix (and the corresponding primal DOF)
%         BC indicates the lagrange DOFs for boundary conditions for each
%            model
%         C1/2 corresponds to the coordinates of the C1/2 matrix (the last
% dimension of the matrix indicates x/y/z)
% each of these index matrix is a N*2 matrix, where the first column
% indicates the first element and the second column indicates the last
% element. The number of lines in K/BC in the number of models, the number
% of lines in C1/2 is the number of models.
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2010
% modif YLG 03/2013: add the mesomicro coupling

% constants
Nm = length(model);
Nc = length(coupling);

% initializations
Nmi = zeros( Nm, 1 );  %taille matrices de raideur
BCi = zeros( Nm, 1 );  %indice debut ddl conditions aux limites -1
Nci = zeros( Nc, 1 );  %taille ddl de couplage
c2m = zeros( Nc, 2 );  %indice numero modeles
x = [];
y = [];
K = [];
z = [];
k = [];
F = [];
indsto = [];

% get sizes and correspondence between coupling and models
for i1 = 1:Nm
    Nmi(i1) = max(model{i1}.K.x);
    BCi(i1) = size(model{i1}.mesh.tri3.X,1);
end
indKi = [ones(Nm,1) BCi];   %indice raideur dans les matrices
indBCi = [indKi(:,end)+1 Nmi];  %indice BC dans les matrice
indK =  indKi+ repmat(cumsum([0;indBCi(1:end-1,2)],1),1,2) ;    %indice raideur dans la matrice globale
indBC =   indBCi+ repmat(cumsum([0;indBCi(1:end-1,2)],1),1,2)   ;     %indice BC dans la matrice globale
%indK = [0 ; cumsum(Nmi(1:Nm)) ];    %
%indBC = indK(1:Nm) + BCi;
for i1 = 1:Nc
    Nci(i1) = max( coupling{i1}.C1.y );
end
for i1 = 1:Nc
    c2m( i1, : ) = coupling{i1}.models;
end
%indBC = [ indBC+1  indK(2:end) ];
%indK = [ indK(1:end-1)+1 indBC(:,1)-1 ];
indC1x = indK(c2m(:,1),:);    %indice ligne couplage matrice globale modele 1
indC = indBC(end,2) + [ [ 1 ; cumsum(Nci(1:end-1)) ] cumsum(Nci) ];   %indice ddl couplage matrice globale
indC2x = indK(c2m(:,2),:);    %indice ligne couplage matrice globale modele 2

% assemble the pure stiffness part
for i1 = 1:Nm
    Ki = model{i1}.K;
    Fi = model{i1}.F;
    x = [ x ; indK(i1,1)-1 + Ki.x ];
    y = [ y ; indK(i1,1)-1 + Ki.y ];
    K = [ K ; Ki.val ];
    z = [ z ; indK(i1,1)-1 + Fi.x ];
    k = [ k ; Fi.y ];
    F = [ F ; Fi.val ];
end


% output
opt = struct( 'K', indK, ...
    'BC', indBC, ...
    'Cy', indC, ...
    'C1x', indC1x, ...
    'C2x', indC2x);

% assemble the coupling parts (including transpose parts)
for i1 = 1:Nc
    C1 = coupling{i1}.C1;
    C2 = coupling{i1}.C2;
    x = [ x ; indC1x(i1,1)-1 + C1.x; indC(i1,1)-1 + C1.y ...
        ; indC2x(i1,1)-1 + C2.x; indC(i1,1)-1 + C2.y ];
    y = [ y ; indC(i1,1)-1 + C1.y; indC1x(i1,1)-1 + C1.x ...
        ; indC(i1,1)-1 + C2.y; indC2x(i1,1)-1 + C2.x];
    K = [ K ; C1.val ; C1.val ; -C2.val ; -C2.val ];
    % warning: this was only checked for a single stochastic model
    
    %%%LAISSER LES CAS ET LES F regarder si ok avec det-sto%
    if or(isfield( model{c2m(i1,1)}.K, 'MC' ),isfield( model{c2m(i1,2)}.K, 'MC' ))
        opt.MC={};
        for ijk=1:length(c2m(i1,:))
            imod1 = c2m( i1, ijk );
            Cmod1 = eval( [ 'coupling{i1}.C' num2str(imod1) ] );
            indCtx(imod1,:) = indK(imod1,:);
            indCt(i1) = indC(i1,end)+1;
            indSc(i1,:) = indC(i1,:);
            indScx(i1) = indC(i1,end)+2;
            xCt=indCtx(imod1,1)-1+Cmod1.xtheta;
            yCt=indCt(i1)*ones(length(xCt),1);
            ySc=indSc(i1,1)-1+Cmod1.xBCpsi;
            xSc=indScx(i1)*ones(length(ySc),1);
            signe = 1*(ijk==1)-1*(ijk==2);
            KCt=signe*Cmod1.Ctheta;
            KSc=Cmod1.BCpsi;
            if ijk==1
                x=[x;xCt;yCt;xSc;ySc];
                y=[y;yCt;xCt;ySc;xSc];
                K=[K;KCt;KCt;KSc;KSc];
            else
                x=[x;xCt;yCt];
                y=[y;yCt;xCt];
                K=[K;KCt;KCt];
            end
            
            if(isfield( model{imod1}.K, 'MC' ))
                indsto(end+1)=imod1;
                Nmc = size(model{imod1}.HomeFE.property,3);
                xKs = model{imod1}.K.x + indK(imod1,1)-1;
                yKs = model{imod1}.K.y + indK(imod1,1)-1;
                Ksi = model{imod1}.K.MC;
                for ijkl = 1:Nmc
                    opt.MC{imod1,ijkl} = sparse(  xKs, yKs, Ksi{ijkl},max(x),max(y));
                end
                
                
            end
        end
    end
end

% create sparse matrix
%K = sparse( x, y, K );
F = sparse( z, k, F, max(x), max(k) );
opt.indsto=indsto;

K=sparse(x,y,K,max(x),max(y));
for ijk=indsto
    K(indK(ijk,1):indBC(ijk,2),indK(ijk,1):indBC(ijk,2))=0;
end
[ xf, yf, F ] = find( F );
F = sparse( xf, yf, F, max(x), max(k) );