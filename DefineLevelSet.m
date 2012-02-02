function LSet = DefineLevelSet( X, weight )
% DEFINELEVELSET to compute the level set function at nodes X
%
% syntax: LSet = DefineLevelSet( mesh, Cint, Cext )
%
%  X1,X2: nodal coordinates matrices
%  weight: structured array containing the fields
%       -'ext': external geometrical limit of the coupling zone  
%               given as a mesh in (X,T) format, with one dimension less  
%               than X1
%       -'int': internal geometrical limit of the coupling zone
%               given in the same format as ext. The curve int should be 
%               located inside the curve ext.
%
%  LSet: structured array containing the fields
%      -'ext': value of the level set function corresponding to the
%              exterior boundary (positive outside the curve, negative
%              inside)
%      -'int': value of the level set function corresponding to the 
%              interior boundary (positive inside the curve, negative 
%              outside)
%
% contact: regis.cottereau@ecp.fr

% R. Cottereau 01/2011

% construct the level set
dint = DistanceMesh2LevelSet( X, weight.int );
dext = DistanceMesh2LevelSet( X, weight.ext );

% change the sign if need be when exterior boundary is fully contained
% inside the interior boundary
if LSinLS( weight.ext, weight.int ) % case 'zoom': coarse model
    dext = -dext;
elseif LSinLS( weight.int, weight.ext ) % case 'zoom': fine model
    dint = -dint;
end

% store level set function
LSet = struct( 'int', dint, ...
               'ext', dext );

%==========================================================================
function dist = DistanceMesh2LevelSet( X, LSet )
% DISTANCEMESH2LEVELSET to compute the distance from a set of points to a
% curve given as a list of points
%
%  syntax: d = DistanceMesh2LevelSet( X, LevelSet )
%
%  X: set of Nn nodes in d dimensions [Nn*d matrix]
%  LevelSet: definition of a curve as a (X,T) pair in (d-1) dimension. The
%            curve must be closed
%
%  dist: list of distances for each node in X [Nn*1 vector] (interior
%        values are negative)

% R. Cottereau 04/2010

% constants
[Nn d] = size( X );

% initializations
dist = zeros( Nn, 1 );

% geometrical dimension
switch d
    
    % 1D case (only two nodes defines the level set, one being +/-Inf)
    case 1
        dist1 = X - LSet.X(1);
        dist2 = X - LSet.X(2);
        dist = min( abs(dist1), abs(dist2) );
        ind = (dist1.*dist2)<0;
        dist(ind) = -dist(ind);

    % 2D case
    case 2
        poli = LSet.X( [LSet.T(:,1); LSet.T(end,2) ], : );
        
        % loop on nodes
        for i1 = 1:Nn
            Xn = X(i1,:);
            dist(i1) = dppoli( Xn, poli );
        end
        
        % change the sign for nodes inside the polygon
        ind = inpolygon( X(:,1), X(:,2), LSet.X(:,1), LSet.X(:,2) );
        dist(ind) = -dist(ind);

    % higher dimensions not implemented yet
    otherwise
        error( 'distance to level set not available in d>2 dimensions' )
end
    
%------------------------------------------------------------------------
function res = dppoli( pto, poli )
% distancia de un punto a un poligono
% pto = [x y]
% poli = solo los nodos. se cierra el primero y el ultimo

n1 = size( poli, 1 ) - 1;
dist1 = zeros(n1,1);
for i = 1:n1
   rec = [ poli(i,1) poli(i,2) ; poli(i+1,1) poli(i+1,2) ];   
   dist1(i) = dpr( pto, rec );
end
res = min( dist1 );


%------------------------------------------------------------------------
function dist = dpr( p, r )
% dist punto recta
% punto = [x y]
% recta = [x1 y1 ; 
%          x2 y2 ]
p1  = r(1,:);
p2  = r(2,:);
dir = p2 - p1;
modd = norma(dir);
if modd ~= 0
   dir = dir / modd;
   pto = p - p1;
   f   = sum( dir.*pto );
   pr  = p1 + dir * f;
   if f < 0 || f > modd
      % proy fuera
      dist = min( dpp( p1, p ), dpp( p2, p ) );
   else
      %proy dentro
      dist = dpp( pr, p );
   end
else
   error( 'MkLS: un vector tiene long cero' )
end

%------------------------------------------------------------------------
function d = dpp( p1, p2 )
% dist punto punto
% punto = [x y]
d = ((p1(1)-p2(1)).^2 + (p1(2)-p2(2)).^2 ).^0.5;

function n = norma( v )
% norma
n = (v(:,1).^2 + v(:,2).^2).^0.5;
