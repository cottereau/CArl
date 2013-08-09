function u = SolveArlequin( K, F, solver, Kmc )
% SOLVEARLEQUIN to solve the coupled system
%
% syntax: sol = SolveArlequin( K, F, coupling )
%
%  K: sparse stiffness matrix
%  F: sparse vector
%  solver: 'direct', or 'FETI' (not implemented yet)
%
%  sol: solution of the system
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% solve depending on the type of solver demanded
switch lower(solver)
    
    % direct solver
    case 'direct'
        if isempty(Kmc)
            u = K \ F;
        else
            Nmc = length( Kmc );
            u = zeros( size(K,1), Nmc );
            for i1= 1:Nmc
                u(:,i1) = (K+Kmc{i1}) \ F;
            end
        end

    % FETI solver
    case 'feti'
        error( 'not implemented yet' )
        
    % unknown solver
    otherwise
        error( 'unknown solver type' )        
end
