function [ Mdl, Cpl ] = SolveArlequin( K, F, Mdl, Cpl, solver, opt )
% SOLVEARLEQUIN to solve the coupled system
%
% syntax: sol = SolveArlequin( K, F, coupling )
%
%  K: sparse stiffness matrix
%  F: sparse vector
%  solver: 'direct', 'montecarlo' or 'FETI' (not implemented yet)
%
%  sol: solution of the system
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2010

% constants
Nm = size(opt.K,1);
Nc = size(opt.Cy,1);

% solve depending on the type of solver demanded
switch lower(solver)
    
    % direct solver
    case 'direct'
        u = K \ F;
        
        % Monte Carlo solver
    case 'montecarlo'
        Nmc = length( opt.MC.Ks );
        % warning: this has only been tested with one stochastic model
        
        % index vectors
        ind0 = opt.K(opt.MC.i,1) : opt.BC(opt.MC.i,2);
        Nmi = length(ind0);
        Nmki = Nmi - diff(opt.BC(opt.MC.i,:)) - 1;
        ind1 = setdiff( opt.K(1,1):opt.BC(end,2), ind0 );
        ind = [ ind1 opt.Cy(1,1):opt.Cy(end,2) opt.Cy(end,2)+2 ];
        ind2 = setdiff( 1:size(K,1), ind );
        
        % submatrices
        K1 = full( K( ind, ind ) );
        K2 = K( ind2, ind2 );
        C = K( ind, ind2 );
        F1 = F( ind );
        F2 = F( ind2 );
        invK1 = pinv( K1 );
        [K1x,K1y,K1val] = find(invK1) ;
        invK1 = sparse(K1x,K1y,K1val) ;
        R = null( full(K1) );
        [Rx,Ry,Rval] = find(R) ;
        R = sparse(Rx,Ry,Rval,size(R,1), size(R,2)) ;
        clear Rx Ry Rval K1x K1y K1val
        f = [F2-C'*invK1*F1; R'*F1];
        c = C'*R;
        
        % MC loop and computations of realizations of the patch
        % solutions
        u2 = zeros( Nmi, Nmc );
        theta = zeros( 1, Nmc );
        alpha = zeros( size(R,2), Nmc );
        for i1= 1:Nmc
            Ks = [ opt.MC.Ks{i1} zeros(Nmi,1); zeros(1,Nmi+1) ];
            k = K2 + Ks - C'*invK1*C;
            u = [k c;c' zeros(size(c,2))]\f;
            u2(:,i1) = u(1:Nmi);
            theta(i1) = u(Nmi+1);
            alpha(:,i1) = u((Nmi+2):end);
        end
        Mdl{opt.MC.i}.uMC = u2(1:Nmki,:);
        Mdl{opt.MC.i}.lambdaBCMC = u2(Nmki+1:Nmi,:);
        Mdl{opt.MC.i}.lambdaThetaMC = theta;
        
        % solution over the deterministic mesh
        u = invK1*(F1-C*mean([u2;theta],2)) + R*mean(alpha,2);
        u = [ u( opt.K(1,1):opt.BC(1,2) );
              mean( u2, 2 )
              u( (length(ind1)+1):end ) ];
        
        % FETI solver
    case 'feti'
        error( 'not implemented yet' )
        
        % unknown solver
    otherwise
        error( 'unknown solver type' )
        
end

% prepare output
for i1 = 1:Nm
    Mdl{i1}.u = full( u( opt.K(i1,1):opt.K(i1,2) ) );
    Mdl{i1}.lambdaBC = full( u( opt.BC(i1,1):opt.BC(i1,2) ) );
end
for i1 = 1:Nc
    Cpl{i1}.lambda = full( u( opt.Cy(i1,1):opt.Cy(i1,2) ) );
end
