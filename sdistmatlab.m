function d = sdistmatlab( X, G, dG1 )
% SDISTMATLAB to compute the distance from a set of nodes X to a closed
% curve G
%
% syntax: d = sdistmatlab( X, G, [dp] )
%
%  X: set of nodes to be checked [Nn*2 matrix]
%  G: piecewise linear curve given by a set of nodes [Ng*2 matrix]. The
%     last node should not be repeated with the first one, althought this
%     curve is implicitely closed
%  dp: distance from all the points in X to all the points in G 
%      [Nn*Ng matrix]. [optional]
%
%  d: distance from all the points in X to each of the segments in G 
%     [Nn*Ng matrix]

% R. Cottereau 06/2010

% constants
Ng = size( G, 1 );
Nx = size( X, 1 );

% minimal distance to nodes of G ( 1 and 2 correspond to the two
% extremities of each segment in G)
if nargin<3
    dG1 = X2X( X, G );
end
dG2 = dG1(:,[2:Ng 1]);

% length, normal vector, and equations for the segments of the curve
lG = x2x( G, G([2:Ng 1],:) );
nG = diff( G([1:Ng 1],:) ) ./ [ lG lG ];
ABC = [ -nG(:,2) nG(:,1) -nG(:,1).*G(:,2)+nG(:,2).*G(:,1) ];

% distance from the nodes in X to each of the segments, considered as an
% infinite line
dl = abs( [ X ones(Nx,1) ] * ABC' );
d = dl;

% choose this last distance only when the node is "between" the two
% extremities of the segment
ind = ( max(dG1,dG2).^2-d.^2 > repmat( lG'.^2, [Nx 1] ) );
d( ind ) = min( dG1(ind), dG2(ind) );

% SUBROUTINES
% distance from a set of nodes to a set of nodes, one by one
function d = x2x( X1, X2 )
d=sqrt( (X2(:,1)-X1(:,1)).^2 + (X2(:,2)-X1(:,2)).^2);
