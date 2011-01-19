function xn=refine1(pts,n)
%REFINE1  Refine grid uniformly around critical points.
%   This is a routine called by FEM1.
%   REFINE1 uniformly refines the grid around points.
%
%   XN = REFINE1(PTS,N),
%   where PTS is a vector containing both endpoints and the points
%   around which the grid will be refined.
%   N is the number of refinements.
%   The number of elements will be
%      #EL = (#PTS-2)*(2*RES-2)
%   where
%      RES = 2^(N-1)+1
%   thus
%      #EL = (#PTS-2)*2^N
%   The number of gridpoints (length(xn)) is then
%      #GRD = #EL + 1
%   For example:
%
%      pts=[-2 0 2]
%
%   yields:
%
%      xn=[-2 0 2]                          for n=0,1
%      xn=[-2 -.6667 0 .6667 2]             for n=2
%      xn=[-2 -1.2 -.6 -.2 0 .2 .6 1.2 2]   for n=3
%
%   See also ADAPT1, FEM1, GENMAT1.

% Copyright (c) 2002-03-13, B. Rasmus Anthin.


error(nargchk(1,6,nargin))
if nargin<2, n=2;end
if n<1, n=1;end
bd=pts([1 end]);
pts=pts(2:end-1);
x=union(vmean(pts),union(pts,bd));
xn=[];
res=2^(n-1)+1;
for i=1:length(x)-1
   xn=union(xn,quadspace(x(i),x(i+1),(-1)^i*res));
end
if n<2,xn=union(bd,pts);end