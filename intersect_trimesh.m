function [ Ti, Xi ] = intersect_trimesh( T1, X1, T2, X2, geps )
% INTERSECT_TRIMESH to create the intersection of two triangular meshes
%
% syntax: [ Ti, Xi ] = intersect_trimesh( T1, X1, T2, X2, geps )
%
%  T,X: connectivity matrices and nodal coordinates for the two input
%       meshes and the output mesh  [Ne*3 matrix] and [Nn*2 matrix]
%  geps: distance below which two points are considered the same [scalar]

% R. Cottereau 06/2010

% limit the area in which to do the research
[ X1, T1, X2, T2 ] = ResizeMesh( X1, T1, X2, T2, geps );

% constants
Ne1 = size(T1,1);
Ne2 = size(T2,1);
Nn1 = size(X1,1);
Nn2 = size(X2,1);

% initialization
Ti = zeros( 2*max(Ne1,Ne2), 3 );
kT = 0;
Xi = zeros( 2*max(Nn1,Nn2), 2 );
kX = 0;

% computation of approximate radii and centers of all the elements
[ c1, r1 ] = CenterRadius( X1, T1 );
[ c2, r2 ] = CenterRadius( X2, T2 );

% computation of distances between each of the two sets of points
d12 = X2X( X1, X2 );

% loop on the elements of the first mesh
for i1=1:Ne1

    Te1 = T1( i1, : );
    Xe1 = X1( Te1, : );
    testint1 = find( test_intersec( c1(i1,:), r1(i1), c2, r2, geps ) );
    
    % loop on the elements of the second mesh
    for i2=1:length( testint1 )
    
        Te2 = T2( testint1(i2), : );
        Xe2 = X2( Te2, : );
%          if i1==4 & testint1(i2)==6
%              keyboard
%          end
        
        % intersection of the two elements
        Xii = intersect_tritri( Xe1, Xe2, geps, 0, d12(Te1,Te2) );
        if isempty( Xii )
            continue
        end
        Nni = size(Xii,1);
        ind = kX+(1:Nni);
        Xi( ind, : ) = Xii;
        kX = kX + Nni;

        % computation of an integration mesh
        Tii = delaunay( Xii(:,1), Xii(:,2), {'Qt','Qbb','Qc','Qz'} );
        Nei = size(Tii,1);
        Ti( kT+(1:Nei), : ) = ind(Tii);
        kT = kT + Nei;

    end
end

% get rid of repeated nodes and renumber the connectivity matrix
Ti = Ti( 1:kT, : );
Xi = Xi( 1:kX, : );
[ Xi, Ti ] = GetRidRepeatedNodes( Xi, Ti );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subfunctions                                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ==== CENTERDIAMETER ==========
function [c,r] = CenterRadius( X, T )
c = ( X(T(:,1),:) + X(T(:,2),:) + X(T(:,3),:) )/2;
r = max( [ (X(T(:,1),1)-c(:,1)).^2+(X(T(:,1),2)-c(:,2)).^2 ...
           (X(T(:,2),1)-c(:,1)).^2+(X(T(:,2),2)-c(:,2)).^2 ...
           (X(T(:,3),1)-c(:,1)).^2+(X(T(:,3),2)-c(:,2)).^2 ], [], 2 );
r = 1.01 * sqrt( r );

% ==== TEST_INTERSEC ==========
function ind=test_intersec( c1, r1, c2, r2, geps )
ind = ( ( (c1(1)-c2(:,1)).^2 + (c1(2)-c2(:,2)).^2 ) < ((r1+r2+geps).^2) );

% ==== REPEATEDNODES ==========
function [X,T] = GetRidRepeatedNodes( X, T )
[X,indI,indJ] = unique( X, 'rows' );
Te = T;
for i1 = length(indJ):-1:1
    T( Te==i1 ) = indJ(i1);
end
% ==== RESIZEMESH ==========
function [ X1, T1, X2, T2 ] = ResizeMesh( X1, T1, X2, T2, geps )
[ xmin1, xmax1 ] = MaxDimMesh( X1, T1 );
[ xmin2, xmax2 ] = MaxDimMesh( X2, T2 );
xmin = max( [ xmin1 ; xmin2 ] ) - geps;
xmax = min( [ xmax1 ; xmax2 ] ) + geps;
[ X1, T1 ] = CutMesh( X1, T1, xmin, xmax );
[ X2, T2 ] = CutMesh( X2, T2, xmin, xmax );
% ==== CUTMESH ==========
function [ X, T ] = CutMesh( X, T, xmin, xmax )
indX = all( [X(:,1)>xmin(1) X(:,2)>xmin(2) ...
             X(:,1)<xmax(1) X(:,2)<xmax(2)], 2 );
indT = any( indX( T ), 2 );
[ X, T ] = ReduceMesh( X, T(indT,:) );
% ==== MAXDIMMESH ==========
function [ xmin, xmax ] = MaxDimMesh( X, T )
d = max( abs( [ X(T(:,1),:)-X(T(:,2),:) ;
                X(T(:,2),:)-X(T(:,3),:) ;
                X(T(:,3),:)-X(T(:,1),:) ] ) );
xmin = min(X)-d;
xmax = max(X)+d;

