% distance from a set of points to a polygon
% dist = dppoli( pto, poli )
function dist = dppoli( pto, poli )
n1 = size( poli, 1 ) - 1;
dist = inf(size(pto,1),1);
for i1 = 1:n1
    dist = min( dist, dpr( pto, poli(i1+(0:1),:) ) );
end
% distance from a set of points to a line
function dist = dpr( p, r )
dist = zeros(size(p,1),1);
p1  = r(1,:);
p2  = r(2,:);
dir = p2 - p1;
modd = norma(dir);
dir = dir / modd;
f = [p(:,1)-p1(1) p(:,2)-p1(2)]*dir';
indneg = f<=0;
dist(indneg) = dpp( p1, p(indneg,:) );
indpos = f>=modd;
dist(indpos) = dpp( p2, p(indpos,:) );
ind = ~indneg & ~indpos;
pint = [p1(1)+dir(1)*f(ind) p1(2)+dir(2)*f(ind)];
dist(ind) = sqrt( sum((pint-p(ind,:)).^2,2) );
% distance between a point p1 and a set of points p2
function d = dpp( p1, p2 )
d = sqrt( (p1(1)-p2(:,1)).^2 + (p1(2)-p2(:,2)).^2 );
% norme of a vector
function n = norma( v )
% norma
n = (v(:,1).^2 + v(:,2).^2).^0.5;

