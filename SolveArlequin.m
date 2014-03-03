function [u,out] = SolveArlequin( model, coupling, solver )
% SOLVEARLEQUIN to solve the coupled system
%
% syntax: sol = SolveArlequin( K, F, coupling )
%
%  K: sparse stiffness matrix
%  F: sparse vector
%  solver: 'direct', or 'latin'
%
%  sol: solution of the system
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% output structure initialization
        out = struct( 'model', {model}, ...
                      'coupling', {coupling} );

% solve depending on the type of solver demanded
switch lower(solver)
    
    % direct solver
    case 'direct'
        [ K, F, opt.iK, opt.iC, Kmc ] = AssembleArlequin( model, coupling );
        if isempty(Kmc)
            u = K \ F;
        else
            Nmc = length( Kmc );
            u = zeros( size(K,1), Nmc );
            for i1= 1:Nmc
                u(:,i1) = (K+Kmc{i1}) \ F;
            end
        end
        out.K = K;
        out.F = F;
        out.opt = opt;

    % LATIN solver
    case 'latin'
        if length(model)>2
            error('latin is implemented only for one coupling');
        end
        K10 = model{1}.K;
        F1 = model{1}.F;
        K20 = model{2}.K;
        F2 = model{2}.F;
        C1 = coupling{1}.C1;
        C1 = [C1; zeros(size(K10,1)-size(C1,1),size(C1,2))];
        C2 = coupling{1}.C2;
        C2 = [C2; zeros(size(K20,1)-size(C2,1),size(C2,2))];
        k1 = 100;
        k2 = 100;
        
        % initialization
        K1 = K10 + k1*(C1*C1');
        U1 = K1\F1;
        phi1 = -k1*C1'*U1;
        K2 = K20 + k2*(C2*C2');
        U2 = K2\F2;
        phi2 = -k2*C2'*U2;
        err = zeros(100,1);
        
        % iterate on U1, U2
        for i1 = 1:1000
            
            % step 1: gluing
            a1 = -phi1 + k1*C1'*U1;
            a2 = -phi2 + k2*C2'*U2;
            w1h = 1/(k1+k2) * (a1+a2);
            w2h = w1h;
            phi1h = k1*w1h - a1;
            phi2h = k2*w2h - a2;
            
            % step 2: uncoupled solutions
            a1h = -phi1h - k1*w1h;
            a2h = -phi2h - k2*w2h;
            U1 = K1\(F1-C1*a1h);
            U2 = K2\(F2-C2*a2h);
            phi1 = -k1*C1'*U1 - a1h;
            phi2 = -k2*C2'*U2 - a2h;
            
            % error indicator
            err(i1) = 2*sqrt(norm(phi1-phi1h)+norm(phi2-phi2h)) ...
                                 / sqrt(norm(phi1+phi1h)+norm(phi2+phi2h));
            
        end
        
        % output
        u = [ U1; U2; phi1h];
        out.opt = struct('iK',[0 size(K10,1) size(K10,1)+size(K20,1)] );
        out.err = err;
    
    % unknown solver
    otherwise
        error( 'unknown solver type' )        
end
