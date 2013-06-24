function [ x, y, C ] = CouplingOperatorFE2D( operator, mesh, opt )

% YLG 06/2013

% constants
[Nne,ne] = size(mesh.Triangulation);

% computation of MassMatrix
m = struct( 'mesh', mesh, ...
            'property', ones(Nne,ne), ...
            'load', zeros(Nne,ne), ...
            'BC', [] );
        elements3 = mesh.Triangulation;
coordinates = mesh.X;
[ M ] = mass_matrixP1_2D_elasticity( elements3,coordinates );
[x,y] = find(M);
indtemp = find(M);
C=M(indtemp);
% choice of the coupling operator
switch operator
    
    % L2 coupling
    case 'L2'            
        
    % H1 coupling
    case 'H1'
lambda=ones([size(elements3,1) 2]); mu=zeros([size(elements3,1) 2]); %Lamme coeficients
%dirichlet = [model.BC.nodes(1:end-1) model.BC.nodes(2:end)];
P = zeros([size(elements3,1) size(elements3,2) 2]);
        [ KK, ~, ~ ] = stifness_matrixP1_2D_elasticity(elements3,coordinates,lambda,mu,P );
        [z,k] = find(KK);
        indtemp = find(KK);
        K = KK(indtemp);
        x = [ z; x ];
        y = [ k; y ];
        C = [ opt.kappa*K; C ];
        
    % unknown coupling operator
    otherwise
        error('unknown coupling operator')
end