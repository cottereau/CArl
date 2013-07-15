function [ Mdl, Cpl ] = SolveArlequin( K, F, Mdl, Cpl, solver, opt )
% SOLVEARLEQUIN to solve the coupled system
%
% syntax: sol = SolveArlequin( K, F, coupling )
%
%  K: sparse stiffness matrix
%  F: sparse vector
%  solver: 'direct', 'montecarlo', 'dmc', or 'FETI' (not implemented yet)
%
%  sol: solution of the system
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% solve depending on the type of solver demanded
switch lower(solver)
    
    % direct solver
    case 'direct'
        u = K \ F;
        
    % monte carlo solver
    case 'montecarlo'
        tic
        Nmc = length( opt.MC.Ks );
        
        % submatrices
        indk = opt.K(1,1) : opt.BC(end,2);
        indc = opt.Cy(1,1) : opt.Cy(end,2)+1;
        C = K( indk, indc );
        K = K( indk, indk );
        nc = size(C,2);
        nk = size(K,1);
        O = spalloc( nc, nc, 0 );
        
        % MC loop and computations of realizations of the patch solutions
        u = zeros( nk+nc, Nmc );
        for i1= 1:Nmc
            u(:,i1) = [K+opt.MC.Ks{i1} C;C' O] \ F;
        end
        toc

    % FETI solver
    case 'feti'
        error( 'not implemented yet' )
        
    % unknown solver
    otherwise
        error( 'unknown solver type' )        
end

% prepare output
for i1 = 1:length(Mdl)
    Mdl{i1}.u = ( u( opt.K(i1,1):opt.K(i1,2), : ) );
    Mdl{i1}.lambdaBC = ( u( opt.BC(i1,1):opt.BC(i1,2), : ) );
end
for i1 = 1:length(Cpl)
    Cpl{i1}.lambda = ( u( opt.Cy(i1,1):opt.Cy(i1,2), : ) );
end