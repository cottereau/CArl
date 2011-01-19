clear all ;
close all ;
clc ;

% load models
load Tests/test1bis.mat ;
model = test ;

% description of the model
X = model.mesh.X ;
T = model.mesh.T ;

% constants
[ Ne nnode ] = size( T ) ;
[ Nn d ] = size( X ) ;

model.alpha = ones(Nn , 1) ;

% property
model.property = model.alpha .* model.property;
E = model.property ;
C = zeros(Nn,1) ;

% vector force
Fb = model.load ;

if d==1
    [A,B,b]=genmat1(X',E',C',Fb') ;
elseif d==2
    [A,B,b]=genmat2(X,T,E,C,Fb) ;
end
Kint = A+B ;
[x,y] = find(Kint) ;
K = full(diag(Kint(x,y))) ;
[z,yz] = find(b) ;
F = full(diag(b(z,yz))) ;
k = ones(size(z,1),1) ;

% add boundary conditions
if ~isempty( model.BC )
    Nbc = length( model.BC.type ) ;
    disp('warning: only dirichlet boundary conditions implemented') ;
    Nx = max(x) ;
    Nf = size(model.load,2) ;
    x = [ x ; Nx+(1:Nbc)'; model.BC.nodes' ] ;
    y = [ y ; model.BC.nodes'; Nx+(1:Nbc)' ] ;
    K = [ K ; ones( 2*Nbc, 1 ) ] ;
    z = [ z; reshape(repmat(Nx+(1:Nbc)',[1 Nf]),Nbc*Nf,1) ] ;
    k = [ k; reshape(repmat(1:Nf,[Nbc 1]),Nbc*Nf,1) ] ;
    F = [ F; reshape(repmat(model.BC.value',[1 Nf]),Nbc*Nf,1) ] ;
end

K = sparse( x, y, K ) ;
F = sparse(z, k, F) ;

% Solution
sol1 = K \ F ;
sol = sol1(1:Nn,1) ;

plot(X,sol,'r')