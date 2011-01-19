function Xi = intersect_tritri( X1, X2, geps, l, d12 )
% INTERSECT_TRITRI: computation of the intersection points between two
% triangular elements
%
%  syntax: Xi = intersect_tritri( X1, X2, geps, lg, [dp] )
%
%  X: nodal coordinates of the two triangles [3*2 matrix]
%     these are implicit triangles , the orientation of which is given by
%     the order in which the nodes are given
%  geps: distance below which two points are considered the same [scalar]
%  dp: distance from all the points in X to all the points in G 
%      [N1*N2 matrix]. [optional]
%
%  Xi: coordinates of the intersection points. The output is not
%      necessarily a triangles so the order of the nodes in Xi is not
%      relevant

% R. Cottereau 06/2010

% distances from the points in X1 to the curve X2
if nargin<5
    d12 = sdistmatlab( X1, X2 );
end
%keyboard
% nodes of X1 inside or exactly on the curve X2
in = inpolygon( X1(:,1), X1(:,2), X2(:,1), X2(:,2) );
on = any( d12 < geps, 2 );
ok = in | on;

% initialiasation
edges = [1 2 ;
         2 3 ;
         3 1];
Xi = X1( ok, : );

% find edges of the triangle X1 that have both nodes ok
type1 = all( ok( edges ), 2 );
if all(type1)
    return
end

% find edges of the triangle X1 that have one node ok and one outside of X2
% and project the outside nodes on the edges of X2
type2 = any( ok( edges ), 2 ) & ~type1;
etype2 = edges( type2, : )';
%itype2 = etype2( ~ok( etype2 ) );
for i1 = 1:size( etype2, 2 )
    XP1 = X1( etype2(:,i1), : );
    [ Xint, lint ] = intersect_2P( XP1, X2(1:2,:) );
    if lint; Xi = [ Xi; Xint ]; end
    [ Xint, lint ] = intersect_2P( XP1, X2(2:3,:) );
    if lint; Xi = [ Xi; Xint ]; end
    [ Xint, lint ] = intersect_2P( XP1, X2([3 1],:) );
    if lint; Xi = [ Xi; Xint ]; end
end

% find nodes of the triangle X1 that have two nodes outside of X2 and 
% perform the double projection on X2
type3 = all( ~ok( edges ), 2 );
etype3 = edges( type3, : )';
for i1 = 1:size( etype3, 2 )
%if ~isempty( etype3 )
    XP1 = X1( etype3(:,i1), : );
    [ Xint, lint ] = intersect_2P( XP1, X2(1:2,:) );
    if lint; Xi = [ Xi; Xint ]; end
    [ Xint, lint ] = intersect_2P( XP1, X2(2:3,:) );
    if lint; Xi = [ Xi; Xint ]; end
    [ Xint, lint ] = intersect_2P( XP1, X2([3 1],:) );
    if lint; Xi = [ Xi; Xint ]; end
end

% find nodes of X2 that are inside the triangle X1 and include them in the
% final list (be careful: NO ORDER !! )
in2 = inpolygon( X2(:,1), X2(:,2), X1(:,1), X1(:,2) );
Xi = [ Xi; X2(in2,:) ];

% find nodes that are repeated and delete them
Xi = unique( Xi, 'rows' );
if size(Xi,1)<3
    Xi = [];
end

%=========================================================================
function [Xi,l] = intersect_2P( X1, X2 )
% intersection of two pairs of points (non-colinear)
v1 = X1(2,:)-X1(1,:);
v2 = X2(2,:)-X2(1,:);
if abs(v1(1)*v2(2)-v1(2)*v2(1))>1e-10
    a = [ -v1' v2' ] \ (X1(1,:)-X2(1,:))';
    Xi = X1(1,:)+a(1)*v1;
    l = all( (a>=0) & a<=1 );
else
    Xi = [];
    l = 0;
end



