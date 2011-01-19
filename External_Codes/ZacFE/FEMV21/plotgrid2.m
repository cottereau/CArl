function plotgrid2(nod2xy,el2nod,fs)
%PLOTGRID  Plots FE nodes and mesh.
%   PLOTGRID2(NOD2XY,EL2NOD[,FONTSIZE])
%   where;
%   NOD2XY=(#NOD)x(XY) = (#NOD)x2
%   EL2NOD=(#EL)x(#NOD)
%   FONTSIZE = 0 supresses numbering.
%   FONTSIZE = 7 by default.

% Copyright (c) 2002-03-31, B. Rasmus Anthin.

if nargin<3, fs=7;end

figure
x=reshape(nod2xy(el2nod,1),size(el2nod));
y=reshape(nod2xy(el2nod,2),size(el2nod));
patch(x',y','y')
if fs
   text(mean(x'),mean(y'),num2str((1:size(el2nod,1))'),'hor','c','fonts',fs)
end

figure
plot(nod2xy(:,1),nod2xy(:,2),'.')
if fs
   text(nod2xy(:,1)',nod2xy(:,2)',num2str((1:size(nod2xy,1))'),'ver','t','fonts',fs)
end
x=nod2xy(:,1);y=nod2xy(:,2);
Dx=max(x)-min(x);
Dy=max(y)-min(y);
axis([min(x) max(x) min(y) max(y)]+[-Dx Dx -Dy Dy]*.1)