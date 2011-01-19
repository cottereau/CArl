function el2nod=genmesh2(nod2xy,boundp,boundm)
%GENMESH2  Generate 2D mesh from grid points.
%   EL2NOD = GENMESH(NOD2XY,BOUND+[,BOUND-]),
%   where EL2NOD is a (#EL)x3 matrix, containing
%   the nodes (cols) for each element (rows), NOD2XY is
%   a (#NOD)x2 matrix containing the coordinates (cols)
%   for each node (rows), BOUND+ is an array containing
%   the nodes along the boundary/hull along some direction
%   (ordered points) and BOUND- is an array (optional) containing
%   the boundary points (ordered) for those points that will be
%   left out from the resulting mesh. With set algebra this
%   can be written as:
%
%   {GRIDPTS} = {(union_i {A_i})/(union_i {B_i}) :
%             : BOUND+_i = D{A_i} , BOUND-_i = D{B_i} , forall i}
%
%   where D is the boundary operator such that D{A_i} and D{B_i} is
%   boundary of set A_i and B_i respectively.
%
%   See also FEM2, PLOTGRID2, DELAUNAY, INPOLYGON.

% Copyright (c) 2002-04-04, B. Rasmus Anthin.
% Revisited 2003-06-16.

error(nargchk(2,3,nargin))
if size(nod2xy,2)~=2, error('NOD2XY must be a (#NOD)x2 matrix.'),end

el2nod=delaunay(nod2xy(:,1),nod2xy(:,2));
T=elmean(nod2xy,el2nod);
Im=[];
Ip=inside(nod2xy,el2nod,boundp);
if nargin==3
   Im=inside(nod2xy,el2nod,boundm);
end
I=setdiff(Ip,Im);
el2nod=el2nod(I,:);

function T=elmean(nod2xy,el2nod)
%make points in the middle of each triangle
T.x=mean(reshape(nod2xy(el2nod,1),size(el2nod)),2);
T.y=mean(reshape(nod2xy(el2nod,2),size(el2nod)),2);

function I=inside(nod2xy,el2nod,bd)
I=[];
bound=bd;
if ~iscell(bound), bound={};bound{1}=bd;end
for i=1:length(bound)
   bound{i}=bound{i}(:)';
   bound{i}=bound{i}([1:end 1]);
   T=elmean(nod2xy,el2nod);
   %check if the triangle midpoints are outside the boundary or not
   isT=inpolygon(T.x,T.y,nod2xy(bound{i},1),nod2xy(bound{i},2));
   I=[I find(isT)'];
end