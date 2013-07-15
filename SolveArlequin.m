function [ Mdl, Cpl ] = SolveArlequin( K, F, Mdl, Cpl, solver, opt, Kmc )
% SOLVEARLEQUIN to solve the coupled system
%
% syntax: sol = SolveArlequin( K, F, coupling )
%
%  K: sparse stiffness matrix
%  F: sparse vector
%  solver: 'direct', 'montecarlo', or 'FETI' (not implemented yet)
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
        keyboard
        N = size(K,1);
        Nmc = length( opt.MC.Ks );
        
        % MC loop and computations of realizations of the patch solutions
        u = zeros( N, Nmc );
        for i1= 1:Nmc
            u(:,i1) = (K+Kmc{i1}) \ F;
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
    Mdl{i1}.u = ( u( (opt.iK(i1)+1):opt.iK(i1+1), : ) );
end
for i1 = 1:length(Cpl)
    Cpl{i1}.lambda = ( u( (opt.iC(i1,3)+1):opt.iC(i1,3), : ) );
end
