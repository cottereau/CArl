function dist = DistanceMesh2LevelSet( X, LevelSet )
% DISTANCEMESH2LEVELSET to compute the distance from a set of points to a
% curve given as a list of points
%
%  syntax: d = DistanceMesh2LevelSet( X, LevelSet )
%
%  X: set of Nn nodes in d dimensions [Nn*d matrix]
%  LevelSet: list of Nl nodes in d dimensions. The last node should be the
%            same as the first one, to close the loop [Nl*d matrix]
%
%  dist: list of distances for each node in X [Nn*1 vector]
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2010

% constants
[Nn d] = size( X );

% initializations
dist = zeros( Nn, 1 );

% geometrical dimension
switch d
    
    % 1D case (only two nodes define the level set)
    case 1
        dist1 = X - LevelSet(1);
        dist2 = X - LevelSet(2);
        dist =  min( abs(dist1), abs(dist2) );
        ind = dist1.*dist2 < 0;
        dist(ind) = - dist(ind);
                
    % 2D case
    case 2
        % loop on nodes
        for i1 = 1:Nn
            Xn = X(i1,:);
            dist(i1) = dppoli( Xn, LevelSet );
        end
        
        % change the sign for nodes inside the polygon
        ind = inpolygon( X(:,1), X(:,2), LevelSet(:,1), LevelSet(:,2) );
        dist(ind) = -dist(ind);

    % higher dimensions not implemented yet
    otherwise
        error( 'distance to level set not available in 3 dimensions' )
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
