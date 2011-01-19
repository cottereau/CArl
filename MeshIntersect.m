function [ Int, Rep] = MeshIntersect( model1, model2, levelSet1, levelSet2 )
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
X1 = model1.mesh.X;
T1 = model1.mesh.T;
X2 = model2.mesh.X;
T2 = model2.mesh.T;

switch model1.mesh.type
    
    % model 1 is of FE type
    case 'FE'
        
        switch model2.mesh.type
            
            % model1 and model2 are both of FE type
            case 'FE'
                % coupling zone for representation purposes
                indT1 = NonConstantAlpha( model1, levelSet1 );
                [ Xr1, Tr1, Xrg1 ] = ReduceMesh( X1, T1(indT1,:) );
                indT2 = NonConstantAlpha( model2, levelSet2 );
                [ Xr2, Tr2, Xrg2 ] = ReduceMesh( X2, T2(indT2,:) );

                % intersect the meshes to get the integration mesh
                disp('warning: model2 is supposed embedded in model1');
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
                M1 = sparse( size(Xi,1), size(X1,1) );
                for i1 = 1:size(Xr1,1)
                    M1( Xr2Xi1{i1}, Xrg1(i1) ) = val1{i1};
                end
                M2 = sparse( size(Xi,1), size(X2,1) );
                for i1 = 1:size(Xr2,1)
                    M2( Xr2Xi2{i1}, Xrg2(i1) ) = val2{i1};
                end
                
                % construction
                
            % model1 is FE and model2 is discrete
            case 'discrete'
                error('not implemented yet')
                
            % unknown type for model 2
            otherwise
                error('unknown type for model 2')
        end
        
    % model 1 is of discrete type
    case 'discrete'
        error('not implemented yet')
        
    % unknown type for model 2
    otherwise
        error('unknown type for model 1');
end

% output
Int = struct( 'X', Xi, 'T', Ti );
Rep{1} = struct( 'type', model1.mesh.type,...
             'X', Xr1, ...
             'T', Tr1, ...
             'M', M1, ...
             'value', {val1}, ...
             'Xr2Xi', {Xr2Xi1}, ...
             'Xrg', Xrg1 );
Rep{2} = struct( 'type', model2.mesh.type,...
             'X', Xr2, ...
             'T', Tr2, ...
             'M', M2, ...
             'value', {val2}, ...
             'Xr2Xi', {Xr2Xi2}, ...
             'Xrg', Xrg2 );


%==========================================================================
function ind = NonConstantAlpha( model, alpha )
% to extract a submesh of model.mesh where the alpha function is not 
% constant
d = size( model.mesh.X, 2 );
T = model.mesh.T;
ind = ( abs(alpha(T(:,1)) - alpha(T(:,2))) > 1e-9 );
if d>1
    ind = ( abs(alpha(T(:,1)) - alpha(T(:,2))) > 1e-9 )...
        | ( abs(alpha(T(:,1)) - alpha(T(:,3))) > 1e-9 )...
        | ( abs(alpha(T(:,2)) - alpha(T(:,3))) > 1e-9 ) ;
end
if d>2
    error('not implemented yet')
end
ind = find( ind );

%==========================================================================

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
disp('warning: this has only been checked in 1D');

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

