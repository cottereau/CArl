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

% R. Cottereau 04/2010
%modif YLG 03/2013: remove condensation for stochastic solver, and add the
%double Monte-Carlo solver

% constants
Nm = size(opt.K,1);
Nc = size(opt.Cy,1);

% solve depending on the type of solver demanded
switch lower(solver)
    
    % direct solver
    case 'direct'
        u = K \ F;
        
        % Monte Carlo solvers
        
    case 'montecarlo2'
        tic
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
        Ftot=[F2;zeros(size(K2,1)-size(F2,1));F1;zeros(size(K1,1)-size(F1,1))];
        invK1 = pinv( K1 );
        [K1x,K1y,K1val] = find(invK1) ;
        invK1 = sparse(K1x,K1y,K1val) ;
        R = null( K1 );
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
            Ksp = spalloc(Nmi,Nmi,Nmi^2) ;
            Ksp(1:Nmki,1:Nmki) = opt.MC.Ks{i1}(1:Nmki,1:Nmki) ;
            Ks = [ Ksp zeros(Nmi,1); zeros(1,Nmi+1) ];
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
        toc
    case 'montecarlo'
        tic
        Nmc = length( opt.MC.Ks );
        % warning: this has only been tested with one stochastic model
        
        % index vectors
        ind0 = opt.K(opt.MC.i,1) : opt.BC(opt.MC.i,2);              %stochastic part index
        Nmi = length(ind0);                                         %length stochastic part
        Nmki = Nmi - diff(opt.BC(opt.MC.i,:)) - 1;                  %length stochastic without BC
        ind1 = setdiff( opt.K(1,1):opt.BC(end,2), ind0 );           %deterministic part index
        ind = [ ind1 opt.Cy(1,1):opt.Cy(end,2) opt.Cy(end,2)+2 ];   %deterministic + coupling deter part index
        ind2 = setdiff( 1:size(K,1), ind );                         %stochastic + coupling part index
        
        % submatrices
        K1 = full( K( ind, ind ) );
        K2 = K( ind2, ind2 );
        C = K( ind, ind2 );
        F1 = F( ind );
        F2 = F( ind2 );
        Ftot=[F2;zeros(size(K2,1)-size(F2,1));F1;zeros(size(K1,1)-size(F1,1))];
        invK1 = pinv( K1 );
        [K1x,K1y,K1val] = find(invK1) ;
        invK1 = sparse(K1x,K1y,K1val) ;
        R = null( K1 );
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
            Ksp = spalloc(Nmi,Nmi,Nmi^2) ;
            Ksp(1:Nmki,1:Nmki) = opt.MC.Ks{i1}(1:Nmki,1:Nmki) ;
            Ks = [ Ksp zeros(Nmi,1); zeros(1,Nmi+1) ];
            k=[K2+Ks C';C K1];
            utot=k\Ftot;
            u2(:,i1)=utot(1:size(K2,1)-1);
            theta(i1) = utot(size(K2,1));
            u1(:,i1)=utot(size(K2,1)+1:end);
        end
        Mdl{opt.MC.i}.uMC = u2(1:Nmki,:);
        Mdl{1}.uMC = u1(opt.K(1,1):opt.BC(1,2)-1,:);
        Mdl{opt.MC.i}.lambdaBCMC = u2(Nmki+1:Nmi,:);
        Mdl{opt.MC.i}.lambdaThetaMC = theta;
        u=mean(u1,2);
        u = [ u( opt.K(1,1):opt.BC(1,2) );
            mean( u2, 2 )
            u( (length(ind1)+1):end ) ];
        toc
        
        
    case 'dmc'
        tic
        %Si plusieurs couplages alors juste une boucle supl?mentaires sur
        %les structures coupling est n?c?ssaire.
        [Nmc imicro] = max([length( opt.MC{1}.Ks ),length(opt.MC{2}.Ks)]);
        [Nmeso imeso] = min([length( opt.MC{1}.Ks ),length(opt.MC{2}.Ks)]);
        
        
        indsto1=opt.K(imeso,1):opt.BC(imeso,2);
        indsto1=[indsto1 opt.Cy(imeso,1):opt.Cy(imeso,2)];
        indsto2=opt.K(imicro,1):opt.BC(imicro,2);
        indu2=opt.K(imicro,1):opt.K(imicro,2);
        indu1=opt.K(imeso,1):opt.K(imeso,2);
        indtheta=opt.BC(imicro,1):opt.BC(imicro,2);
        indlambda=opt.BC(imeso,1):opt.BC(imeso,2);
        indPsi=opt.Cy(1,1):opt.Cy(1,2);
        
        C  = full(K(indsto1,indsto2));
        K1 = full(K(indsto1,indsto1));
        K2 = full(K(indsto2,indsto2));
        F1 = F( indsto1 );
        F2 = F( indsto2 );
        Ftot=[F2;F1];
        theta=[];
        for i1= 1:Nmc
            Ksp2(1:length(indsto2)-length(indtheta),1:length(indsto2)-length(indtheta)) = opt.MC{imicro}.Ks{i1}(1:length(indsto2)-length(indtheta),1:length(indsto2)-length(indtheta)) ;
            Ks2 = [ Ksp2 zeros(length(indu2),length(indtheta)); zeros(length(indtheta),length(indsto2)) ];
            Ksp1(1:length(indu1),1:length(indu1)) = opt.MC{imeso}.Ks{ceil(i1/(Nmc/Nmeso))}(1:length(indu1),1:length(indu1)) ;
            Ks1 = [ Ksp1 zeros(length(indu1),length(indlambda)+length(indPsi)); zeros(length(indlambda)+length(indPsi),length(indsto1)) ];
            k=[K2+Ks2 C';C K1+Ks1];
            utot=k\Ftot;
            u2(:,i1)=utot(1:length(indu2));
            theta = [theta utot(length(indu2)+1:length(indu2)+length(indtheta))];
            u1(:,i1)=utot(length(indsto2)+1:end);
        end
        Mdl{opt.MC{2}.i}.uMC = u2;
        Mdl{opt.MC{1}.i}.uMC = u1(1:length(indu1),:);
        Mdl{opt.MC{2}.i}.lambdaBCMC = u1(length(indu1)+1:end,:);
        Mdl{opt.MC{2}.i}.lambdaThetaMC = theta;
%         u=mean(u1,2);
%         u = [ u( 1:length() );
%             mean( u2, 2 )
%             u( (length(indu1)+length(indlambda)+1):end ) ];
        toc
        % FETI solver
    case 'feti'
        error( 'not implemented yet' )
        
        % unknown solver
    otherwise
        error( 'unknown solver type' )
        
end

% prepare output
% for i1 = 1:Nm
%     Mdl{i1}.u = full( u( opt.K(i1,1):opt.K(i1,2) ) );
%     Mdl{i1}.lambdaBC = full( u( opt.BC(i1,1):opt.BC(i1,2) ) );
% end
% for i1 = 1:Nc
%     Cpl{i1}.lambda = full( u( opt.Cy(i1,1):opt.Cy(i1,2) ) );
% end