function [x,y,K,z,F] = K_2D_elasticity(model)
% DOF = [U11 U21 U12 U22 U13 U23 ... U1N U2N]
%addpath('.\library_vectorization')     %in windows
%addpath('./library_vectorization/')     %in linux

%levels=0; %maximum uniform refinement level
E = model.property;
nu = model.HomeFE.poisson;
lambda=E*nu./((1+nu)*(1-2*nu)); mu=E/(2*(1+nu)); %Lamme coeficients

elements3 = model.mesh.Triangulation;
coordinates = model.mesh.X;
%dirichlet = [model.BC.nodes(1:end-1) model.BC.nodes(2:end)];
P = model.load;


    %[K,areas]=stifness_matrixP1_2D_elasticity(elements3,coordinates);
    [stif areas,force]=stifness_matrixP1_2D_elasticity(elements3,coordinates,lambda,mu,P);
[x,y] = find(stif);
indtmp = find(stif);
K = full(stif(indtmp));
z=find(force);
F=force(z); %IMPLEMENTATION A VERIFIER


    
if isfield(model.HomeFE,'BC')&&~isempty( model.HomeFE.BC )
    
    % Dirichlet Boundary Conditions
    ind = find( model.HomeFE.BC.type == 'U' );
    Nbc = length( ind );
    Nx = max(x);
    Nf = 1;
    x = [ x ; Nx+(1:2*Nbc)'; reshape([2*model.HomeFE.BC.nodes(ind)'-1;2*model.HomeFE.BC.nodes(ind)'],length(ind)*2,1) ];
    y = [ y ; reshape([2*model.HomeFE.BC.nodes(ind)'-1;2*model.HomeFE.BC.nodes(ind)'],length(ind)*2,1); Nx+(1:2*Nbc)' ];
    K = [ K ; ones( 4*Nbc, 1 ) ];
    z = [ z; reshape(Nx+(1:2*Nbc)',2*Nbc*Nf,1) ];
    F = [ F; reshape(model.HomeFE.BC.value(ind,:)',2*Nbc*Nf,1) ];
    
    % Neumann Boundary Conditions
    ind = find( model.HomeFE.BC.type == 'F' );
    Nbc = length( ind );
    z = [ z; reshape([2*model.HomeFE.BC.nodes(ind)'-1;2*model.HomeFE.BC.nodes(ind)'],length(ind)*2,1) ];
    F = [ F; reshape(repmat(model.HomeFE.BC.value(ind)',[1 Nf]),2*Nbc*Nf,1) ];
    
end