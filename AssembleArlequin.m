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

% constants
Nm = length(model);
Nc = length(coupling);

% initializations
Nmi = zeros( Nm, 1 );
BCi = zeros( Nm, 1 );
Nci = zeros( Nc, 1 );
c2m = zeros( Nc, 2 );
x = [];
y = [];
K = [];
z = [];
k = [];
F = [];

% get sizes and correspondence between coupling and models
for i1 = 1:Nm
    Nmi(i1) = max(model{i1}.K.x);    
    BCi(i1) = size(model{i1}.mesh.X,1);
end
indK = [0 ; cumsum(Nmi(1:Nm)) ];
indBC = indK(1:Nm) + BCi;
for i1 = 1:Nc
    Nci(i1) = max( coupling{i1}.C1.y );
end
for i1 = 1:Nc
    c2m( i1, : ) = coupling{i1}.models;
end
indBC = [ indBC+1  indK(2:end) ];
indK = [ indK(1:end-1)+1 indBC(:,1)-1 ];
indC1x = indK(c2m(:,1),:);
indC = indBC(end,2) + [ [ 1 ; cumsum(Nci(1:end-1)) ] cumsum(Nci) ];
indC2x = indK(c2m(:,2),:);

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

% assemble the coupling parts (including transpose parts)
for i1 = 1:Nc
    C1 = coupling{i1}.C1;
    C2 = coupling{i1}.C2;
    x = [ x ; indC1x(i1,1)-1 + C1.x; indC(i1,1)-1 + C1.y ...
            ; indC2x(i1,1)-1 + C2.x; indC(i1,1)-1 + C2.y ];
    y = [ y ; indC(i1,1)-1 + C1.y; indC1x(i1,1)-1 + C1.x ...
            ; indC(i1,1)-1 + C2.y; indC2x(i1,1)-1 + C2.x];
    K = [ K ; C1.val ; C1.val ; -C2.val ; -C2.val ];
end

% create sparse matrix
K = sparse( x, y, K );
F = sparse( z, k, F, max(x), max(k) );

% output
opt = struct( 'K', indK, ...
              'BC', indBC, ...
              'Cy', indC, ...
              'C1x', indC1x, ...
              'C2x', indC2x  );

% Monte Carlo case
for i1 = 1:Nc
    if strcmp( coupling{i1}.mediator.type, 'stochastic' )
        
        % warning: this was only checked for a single stochastic model
        
        % constants
        isto = c2m( i1, 1+isfield( model{c2m(i1,2)}.K, 'MC' ) );
        Csto = eval( [ 'coupling{i1}.C' num2str(isto) ] );
        indCstox = eval( [ 'indC' num2str(isto) 'x' ] );
        idet = setdiff( c2m(i1,:), isto );
        Cdet = eval( [ 'coupling{i1}.C' num2str(idet) ] );
        indCdetx = eval( [ 'indC' num2str(idet) 'x' ] );


        % adding the additional constraints in the global stiffness matrix
        [ x0, y0, K0 ] = find( K );
        n = max(x0);
        inddet = indCdetx(1) - 1 + Cdet.xtheta ;
        indsto = indCstox(1) - 1 + Csto.xtheta ;
        indy = indC(i1,1) - 1 + Cdet.xBCpsi;
        ndet = (1+n) * ones( size(inddet) );
        nsto = (1+n) * ones( size(indsto) );
        ny = (2+n) * ones( size(indy) );
        x0 = [ x0; inddet; ndet; indsto; nsto; indy; ny ];
        y0 = [ y0; ndet; inddet; nsto; indsto; ny; indy ];
        K0 = [ K0; Cdet.Ctheta; Cdet.Ctheta;
                   -Csto.Ctheta; -Csto.Ctheta; 
                   Cdet.BCpsi; Cdet.BCpsi ];
        K = sparse( x0, y0, K0 );

        % increase the size of the force matrix
        [ x, y, F ] = find( F );
        F = sparse( x, y, F, size(K,1), 1 );
        
        % erase the contribution of the stochastic model
        ind = indK(isto,1) : indK(isto,2);
        K( ind, ind ) = 0;

        % construction of stochastic matrices with indices
        Ki = model{isto}.K;
        Nmc = length( Ki.MC );
        Ksi = cell( Nmc, 1 );
        xKs = indK(isto,1) : indK(isto,2);
        for i2 = 1:Nmc
            Ksi{i2} = sparse( Ki.x, Ki.y, Ki.MC{i2} );
        end
        opt.MC = struct( 'i', isto, 'K0', K0, 'xKs', xKs, 'Ks', {Ksi} );

    end
end
