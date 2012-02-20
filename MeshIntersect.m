function [ Int, Rep] = MeshIntersect( mesh1, mesh2, LSet1, LSet2 )
% MESHINTERSECT to create the mesh at the intersection between two meshes,
% both in terms of support (for the definition of the interpolation
% functions) and for integration purposes
%
%  syntax: [ci,c1,c2] = MeshIntersect( model1, model2 )
%
%  model1, model2 : structured arrays containing the information relative
%         to the models, in particular:
%       - 'type': describes the type of representation used for the
%                 geometry. This is used for all interpolation purposes,
%                 in particular for the definition of weight functions
%                 Implemented: {'FE' 'discrete'}
%       - 'mesh': array dependent on 'type'.
%
%  ci: mesh used for integration purposes. It is a structured array
%     containing the fields 'X' and 'T'
%  c1,c2: meshes used for representation purposes for each of the two
%         models, and passage between representation and integration
%         meshes. Structured arrays containing the fields
%       - 'Tr' describes the mesh for representation purposes by giving the
%              appropriate indices of model.mesh.T
%       - 'r2i' cell giving the list of elements of the integration mesh Ti
%               that are in each element of the representation mesh Tr
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2010

% constants
X1 = mesh1.X;
T1 = mesh1.Triangulation;
X2 = mesh2.X;
T2 = mesh2.Triangulation;
d = size( X1, 2 );

% coupling zone for representation purposes
indT1 = DefineCouplingElements( T1, LSet1.int.*LSet1.ext );
[ Xr1, Tr1, Xrg1 ] = ReduceMesh( X1, T1(indT1,:) );
indT2 = DefineCouplingElements( T2, LSet2.int.*LSet2.ext );
[ Xr2, Tr2, Xrg2 ] = ReduceMesh( X2, T2(indT2,:) );

% intersect the meshes to get the integration mesh
% warning: model2 is supposed embedded in model1;
Xi = X2;
Ti = T2(indT2,:);
[ Xi, Ti ] = ReduceMesh( Xi, Ti );

% compute the passage matrices in terms of elements
Tr2Ti1 = TR2TI( Tr1, Xr1, Xi, Ti );
Tr2Ti2 = TR2TI( Tr2, Xr2, Xi, Ti );

% compute the passage matrices in terms of nodes
% = get the values of the basis functions in the
% representation meshes at the nodes of the integration
% mesh
[Xr2Xi1,val1] = XR2XI( Xr1, Tr1, Xi, Ti, Tr2Ti1 );
[Xr2Xi2,val2] = XR2XI( Xr2, Tr2, Xi, Ti, Tr2Ti2 );

% construction of passage matrices from the integration
% meshes to the representation meshes
Ni = size(Xi,1);
m1x = zeros(Ni,1); m1y = zeros(Ni,1); m1v = zeros(Ni,1);
for i1 = 1:size(Xr1,1)
    Ni1 = length(Xr2Xi1{i1});
    ind = (i1-1)*Ni + (1:Ni1);
    m1x(ind) = Xr2Xi1{i1};
    m1y(ind) = ones(Ni1,1)*Xrg1(i1);
    m1v(ind) = val1{i1};
end
ind = (m1x~=0);
M1 = sparse( m1x(ind), m1y(ind), m1v(ind), Ni, size(X1,1) );
m2x = zeros(Ni,1); m2y = zeros(Ni,1); m2v = zeros(Ni,1);
for i1 = 1:size(Xr2,1)
    Ni2 = length(Xr2Xi2{i1});
    ind = (i1-1)*Ni + (1:Ni2);
    m2x(ind) = Xr2Xi2{i1};
    m2y(ind) = ones(Ni2,1)*Xrg2(i1);
    m2v(ind) = val2{i1};
end
ind = (m2x~=0);
M2 = sparse( m2x(ind), m2y(ind), m2v(ind), Ni, size(X2,1) );
% M1 = sparse( size(Xi,1), size(X1,1) );
% for i1 = 1:size(Xr1,1)
%     M1( Xr2Xi1{i1}, Xrg1(i1) ) = val1{i1};
% end
% M2 = sparse( size(Xi,1), size(X2,1) );
% for i1 = 1:size(Xr2,1)
%     M2( Xr2Xi2{i1}, Xrg2(i1) ) = val2{i1};
% end            
        
% output
if d==1
    Int.mesh = struct( 'Triangulation', Ti, 'X', Xi );
    mesh1 = struct( 'Triangulation', Tr1, 'X', Xr1 );
    mesh2 = struct( 'Triangulation', Tr2, 'X', Xr2 );
elseif d==2
    Int.mesh = TriRep( Ti, Xi(:,1), Xi(:,2) );
    mesh1 = TriRep( Tr1, Xr1(:,1), Xr1(:,2) );
    mesh2 = TriRep( Tr2, Xr2(:,1), Xr2(:,2) );
end
Rep{1} = struct( 'mesh', mesh1, ...
                 'M', M1, ...
                 'value', {val1}, ...
                 'Xr2Xi', {Xr2Xi1}, ...
                 'Xrg', Xrg1 );
Rep{2} = struct( 'mesh', mesh2, ...
                 'M', M2, ...
                 'value', {val2}, ...
                 'Xr2Xi', {Xr2Xi2}, ...
                 'Xrg', Xrg2 );

%==========================================================================
function ind = DefineCouplingElements( T, LSp )
% to extract a submesh of model.mesh where the product of level set
% functions is negative (ie elements in the coupling zone)
ind0 = abs(LSp) <= 1e-9;
ind = find( all( (LSp(T)>=1e-9) | ind0(T), 2 ) );

%==========================================================================
function indT = TR2TI( Tr, Xr, Xi, Ti )
% for each element in Tr, finds the elements of Ti that are inside

% constants
Ne = size(Tr,1);
d = size(Xr,2);

% initialization
indT = cell(Ne,1);

% loop on representation elements
for i1 = 1:Ne

    % 1D case
    if d==1
        ind = find( (( Xi - Xr(Tr(i1,1)) > -eps ) ) ...
                  & (( Xi - Xr(Tr(i1,2)) <  eps ) ) );
        
    % 2D case
    elseif d==2
        Xv = Xr( Tr(i1,:), 1 );
        Yv = Xr( Tr(i1,:), 2 );
        % parameter 1e-6
        param = 1e-6;
        in1 = inpolygon(Xi(:,1),Xi(:,2),Xv+param,Yv) ;
        in2 = inpolygon(Xi(:,1),Xi(:,2),Xv-param,Yv) ;
        in3 = inpolygon(Xi(:,1),Xi(:,2),Xv,Yv+param) ;
        in4 = inpolygon(Xi(:,1),Xi(:,2),Xv,Yv-param) ;
        in5 = inpolygon(Xi(:,1),Xi(:,2),Xv,Yv) ;
        in = in1 | in2 | in3 | in4 | in5 ;  
        ind = find(in) ;
        
    % not implemented for 3D
    else
        error('not implemented for 3D')
    end

    indT{i1} = find( all( ismember(Ti,ind), 2 ) );
    
end

%==========================================================================
function [Xr2Xi,val] = XR2XI( Xr, Tr, Xi, Ti, Tr2Ti )
% for each node in Xr, find the nodes in Xi that are inside the elements
% that touch it, and the value of the linear FE basis function centered on
% Xr
% warning: this has only been checked in 1D

% constants
[ Nnr d ] = size(Xr); 

% initialization
val = cell( Nnr, 1 );
Xr2Xi = cell( Nnr, 1 );

% loop on the nodes in Xr
for i1 = 1:Nnr
    
    % find the elements in contact with this node
    [ind,indj] = find( Tr == i1 );
    
    % constants and initialization
    Nind = length( ind );
    vali = cell(Nind,1);
    Xr2Xii = cell(Nind,1);
    
    % loop on these elements to get the corresponding nodes and their
    % position in the unit element
    for i2 = 1:Nind
        tmp = unique( Ti( Tr2Ti{ind(i2)}, : ) );
        Xr2Xii{i2} = tmp(:);
        Xloci = Xi( Xr2Xii{i2}, : );
        Xrloc = Xr(Tr(ind(i2),:),:);
        Xloc = GetLocalCoordinates( Xloci, Xrloc );
        if d == 1
            N = [ 1-Xloc Xloc ];
        elseif d == 2
            N = [ 1-Xloc(:,1)-Xloc(:,2) Xloc ];
        else
            error( 'not implemented yet' )
        end
        vali{i2} = N(:,indj(i2));
    end
    
    % get rid of repeated nodes
    [ Xr2Xi{i1}, indk ] = unique( cat( 1, Xr2Xii{:} ) );
    vali = cat( 1, vali{:} );
    val{i1} = vali( indk );
    
end

%==========================================================================
function Xl = GetLocalCoordinates( X, Xe )
% to compute the coordinates of X in the local reference element given by
% Xe;

% constants
d = size(X,2);

if d==1
    Xl = (X-Xe(1))/(Xe(2)-Xe(1));
elseif d==2
    X = [ X(:,1)-Xe(1,1) X(:,2)-Xe(1,2) ];
    P = [-1 1 0; -1 0 1] * Xe;
    Xl = X/P;
end

