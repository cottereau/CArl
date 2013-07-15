function [ K, F, opt, Kmc ] = AssembleArlequin( model, coupling )
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

% constants
Nm = length(model);    % number of models
Nc = length(coupling); % number of coupling pairs

% size of the block matrices and global matrix
Nki = zeros( Nm, 1 );  % size of stiffness matrices (including BC)
Nci = zeros( Nc, 1 );  % size of coupling matrices
c2m = zeros( Nc, 2 );  % numbering of coupling pairs
for i1 = 1:Nm
    Nki(i1) = size(model{i1}.K,1);
end
for i1 = 1:Nc
    Nci(i1) = size( coupling{i1}.C1, 2 );
    c2m(i1,:) = coupling{i1}.models;
end
N = sum(Nki) + sum(Nci);

% correspondance between
% indices of the block matrices in the global matrix
Nci = sum(Nki) + [0;cumsum(Nci)];
iK = [0;cumsum(Nki)]; % first indices within the global matrix
iC = zeros(Nc+1,3);   % first indices within the globla matrix [ix1 ix2 iy]
for i1 = 1:Nc
    iC(i1,:) = [ iK(coupling{i1}.models)' Nci(i1) ];
end
iC(end,end) = N;

% initialization
K = sparse( N, N );
F = sparse( N, 1 );

% assemble the global matrix from the block matrices
for i1 = 1:Nm
    [ x, y, Ki ] = find( model{i1}.K );
    K = K + sparse( iK(i1)+x, iK(i1)+y, Ki, N, N );
    [ x, ~, Fi ] = find( model{i1}.F );
    F = F + sparse( iK(i1)+x, 1, Fi, N, 1 );
end
for i1 = 1:Nc
    [ x, y, C1 ] = find( coupling{i1}.C1 );
    C1 = sparse( iC(i1,1)+x, iC(i1,3)+y, C1, N, N );
    [ x, y, C2 ] = find( coupling{i1}.C2 );
    C2 = sparse( iC(i1,2)+x, iC(i1,3)+y, C2, N, N );
    K = K + C1 + C1' - C2 - C2';
end

% constructing the MonteCarlo block matrices when needed
Kmc = [];
for i1 = 1:Nm
    if ~isempty( model{i1}.Kmc );
        Nmc = length(model{i1}.Kmc);
        if ~exist( 'Kmc', 'var' )
            Kmc = cell(Nmc,1);
            for i2 = 1:Nmc
                Kmc{i2} = sparse( N, N );
            end
        end
        for i2 = 1:Nmc
            [ x, y, Ki ] = find( model{i1}.Kmc{i2} );
            Kmc{i2} = Kmc{i2} + sparse( iK(i1)+x, iK(i1)+y, Ki, N, N );
        end
    end
end

% output
opt = struct( 'iK', iK, 'iC', iC );
