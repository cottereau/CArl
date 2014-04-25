function [u,out] = SolveArlequin( model, coupling, solver, out )
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
        
        % parameters
        errmax = 1e-2;
        nitermax = 1e3;
        relax = .8;
        
        % operators
        K1 = model{1}.K;
        F1 = model{1}.F;
        K2 = model{2}.K;
        F2 = model{2}.F;
        C1 = coupling{1}.C1;
        C1 = [C1; zeros(size(K1,1)-size(C1,1),size(C1,2))];
        C2 = coupling{1}.C2;
        C2 = [C2; zeros(size(K2,1)-size(C2,1),size(C2,2))];
        M = coupling{1}.M;
        P1 = M\C1';
        P2 = M\C2';
        
        % search directions
        lambda = model{2}.FE2D.lambda;
        mu = model{2}.FE2D.mu;
        young = mu*(3*lambda+2*mu)/(lambda+mu);
        gluelength = -2*min(coupling{1}.domain.dist(coupling{1}.domain.mesh.Points));
        k =  young/gluelength;
        %k = 160;
        
        ku1 = k;
        ku2 = k;
        kd1 = k;
        kd2 = k;
        H1 = K1 + kd1*(C1*(M\C1'));
        H2 = K2 + kd2*(C2*(M\C2'));
        
        % initialization
        U1 = H1\F1;
        U2 = H2\F2;
        w1 = P1*U1;
        w2 = P2*U2;
        phi1 = -kd1*w1;
        phi2 = -kd2*w1;
        
        indic = [ 1 zeros(1,nitermax-1) ];
        i1 = 1;     
        
        while (indic(i1) > errmax)&&(i1 <= nitermax-1)
            
            i1 = i1 + 1;
            
            % step 1: gluing
            a1 = -phi1 + ku1*w1;
            a2 = -phi2 + ku2*w2;
            w1h = 1/(ku1+ku2) * (a1+a2);
            w2h = w1h;
            phi1h = ku1*w1h - a1;
            phi2h = ku2*w2h - a2;
            
            % step 2: uncoupled solutions
            a1h = -phi1h - kd1*w1h;
            a2h = -phi2h - kd2*w2h;
            U1old = U1;
            U2old = U2;
            U1 = H1\(F1-C1*a1h);
            U2 = H2\(F2-C2*a2h);
            w1 = P1*U1;
            w2 = P2*U2;
            phi1 = -kd1*w1 - a1h;
            phi2 = -kd2*w1 - a2h;
            
            % relaxation
            U1 = relax*U1+(1-relax)*U1old;
            U2 = relax*U2+(1-relax)*U2old;
            
            % error indicator
            indic(i1) = 2*sqrt(norm(phi1-phi1h)^2+norm(phi2-phi2h)^2) ...
                / sqrt(norm(phi1+phi1h)^2+norm(phi2+phi2h)^2);
                        
        end
                
        % output
        u = [ U1; U2; phi1h ];
        out.opt = struct('iK',[0 size(K1,1) size(K1,1)+size(K2,1)] );
        out.indic = indic;
        
        % unknown solver
    otherwise
        error( 'unknown solver type' )
end
out.U = u;
