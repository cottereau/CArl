function [xnh,ah,bh,sh]=adapt1(xn0,tol,alpha,beta,source,bd)
%ADAPT1  Generate adaptively refined grid.
%   [XN,ALPHA,BETA,S] = ADAPT1(XN0,TOL,ALPHA,BETA,S,BD),
%   where XN is the adaptively generated grid from XN0, ALPHA, BETA and S
%   (see FEM1 for more info).
%   N is the number of refinements and PTS is a vector containing
%   the points around which refinement is made. See REFINE1 for more
%   info on the number of elements returned and how the parameter
%   N affects the resulting resolution of the grid.
%
%   See also REFINE1, FEM1, GENMAT1.

% Copyright (c) 2002-03-24, B. Rasmus Anthin.

a=ev(alpha,xn0);
b=ev(beta,xn0);
s=ev(source,xn0);
U=fem1(xn0,a,b,s,bd).';

xn1=union(vmean(xn0),xn0);
ah=ev(alpha,xn1);
bh=ev(beta,xn1);
sh=ev(source,xn1);
Uh=fem1(xn1,ah,bh,sh,bd).';

xnh=xn0(1);
dx0=diff(xn0(1:2));
i0=1:2:length(xn1);
while 1
   xi1=xnh(end);
   if xi1<xn0(1),error('x_i+1<x0.'),end
   I=vfind(xn1,xi1);
   xs=stagger(xn1);
   Uhs=stagger(Uh);
   dUh=diff(Uhs)./diff(xs);
   ahs=stagger(ah);
   daUh=diff(ahs.*Uhs)./diff(xs);
   S=int(xn0,s);
   Sh=int(xn1,sh);
   xi=(tol-L2norm(xn0,S-Sh(i0),1./a)-L2norm(xn0,(a-ah(i0)).*dUh(i0),1./a)-L2norm(xn0,(b-bh(i0)).*Uh(i0),1./a))/...
      max(daUh(I:end)+sh(I:end))+...
      xi1;
   if xi>xn0(end), break;end
   xnh=[xnh xi];
end
if xnh(end)<xn0(end),xnh=[xnh xn0(end)];end
ah=ev(alpha,xnh);
bh=ev(beta,xnh);
sh=ev(source,xnh);

%---------------------------------------------

function y=stagger(x)
mid=vmean(x);
bd=2*x([1 end])-mid([1 end]);
y=[bd(1) mid bd(2)];

function N=Enorm(x,f,a)
xs=stagger(x);
fs=stagger(f);
N=sqrt(trapz(x,a.*(diff(fs)./diff(xs)).^2));

function N=L2norm(x,f,a)
N=sqrt(trapz(x,a.*f.^2));

function I=int(x,f)
for i=2:length(x)
   I(i)=trapz(x(1:i),f(1:i));
end

function y=ev(f,x)
if isa(f,'inline') | exist(f)==2
   y=f(x);
else
   y=eval(f);
end
if length(y)==1, y=y*ones(size(x));end

function y=fixends(x)
y=x;
y([1 end])=2*x([2 end-1])-x([3 end-2]);

function i=vfind(x,xv)
x2=vmean(x);
x2=[x(1) x2 x(end)];
i=find(x2(1:end-1)<=xv & xv<x2(2:end));
if xv==x(end),i=length(x);end