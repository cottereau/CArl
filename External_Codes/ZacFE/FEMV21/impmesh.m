function [nod2xy,el2nod,bound,W,G]=impmesh(p,e,t)
%[NOD2XY,EL2NOD,BOUND,W,G] = IMPMESH(P,E,T)
%   P E T are variables exported from PDETOOL
%   NOD2XY and EL2NOD are used by FEM2
%   BOUND contains the boundary nodes
%   W contains all of the points in the set
%   G contains only the boundary points
%
%   See also PDETOOL.

% Copyright (c) 2002-04-14, B. Rasmus Anthin.

nod2xy=p';
el2nod=t(1:3,:)';
e=e(1:2,:);
bound=reshape(e,[1 prod(size(e))])';
bound=bound(1:end-1);
W.x=nod2xy(:,1);W.y=nod2xy(:,2);
G.x=nod2xy(bound,1);G.y=nod2xy(bound,2);
