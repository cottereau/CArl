%2D mesh
%close all
clear all
clc
couplage='zoom';
Nmcmeso=5;
Nmcmicro=5;
%Nx=31;
%Ny=15;

%maillage grossier

A1=-1.;
A2=1.;
B1=0.;
B2=1.;
lx=0.1;
ly=0.1;
a=0.5;%profondeur indentation
b=0.5;%demi largeur indentation

%maillage fin
ap=a;
bp=b;
am=0.5;
bm=1.25*lx;
A1m=-bm;
A2m=bm;
B1m=am-1.25*ly;
B2m=am+1.25*ly;
nm=1;

L = 0.05;
dL=0.5;

%lxfin=lx/4/4;
%lyfin=ly/4/4;

%rafinnement
A1r = A1m-lx/4;%-0.1-lx/2;
A2r = A2m+lx/4;%0.1+lx/2;
B1r = B1m-lx/4;%0.4-lx/2;
B2r = B2m+lx/4;%0.6+lx/2;

%couplage

cx1=-2*lx/4;
cx2=2*lx/4;
cy1=ap-ly/2;
cy2=ap+ly/2;

% %test mesh2d
% 
% %ligne 0
% 
% l0x = (cx2:lx/2:b)';
% l0y = (cy2:ly/2:B2)';
% l0 = [l0x l0y];
% 
% %ligne 1
% l1x = (b:lx:A2)';
% l1=[l1x B2*ones(size(l1x))];
% 
% %ligne 2
% l2y = (B2:-ly:B1)';
% l2 = [A2*ones(size(l2y)) l2y];
% 
% %ligne 3
% l3x = (A2:-lx:0)';
% l3 = [l3x B1*ones(size(l3x))];
% 
% %ligne 4
% l4y = (0:ly:B1m)';
% l4=[0*ones(size(l4y)) l4y];
% 
% %ligne 5
% l5x = (0:lx:A2m)';
% l5=[l5x B1m*ones(size(l5x))];
% 
% %ligne 6
% l6y = (B1m:ly:B2m)';
% l6 = [A2m*ones(size(l6y)) l6y];
% 
% % %ligne 7
% % l7x = (A2r:-lx:b)';
% % l7 = [l7x B2r*ones(size(l7x))];
% 
% %ligne 8 
% 
% l8x = [A2m:-lx:0]';
% l8y = [B2m:-ly:a]';
% l8 = [l8x l8y];
% 
% %ligne 9 
% 
% l9y = (a:-ly:B1m)';
% l9=[0*ones(size(l9y)) l9y];
% 
% l=[l0;l1;l2;l3;l4;l5;l6;l8;l9];
% figure;plot(l(:,1),l(:,2))
% [p,t,stats] = mesh2d(l);
% TRItest = TriRep(t,p);
% triplot(TRItest)

%maillage meso

Xgr=A1:lx:A2;%linspace(A1,A2,Nx);
Ygr=B1:ly:B2;%linspace(B1,B2,Ny);
%Xgr2=-0.2:lx/2:0.2;%linspace(A1,A2,Nx);
%Ygr2=0.3:ly/2:0.7;%linspace(B1,B2,Ny);
Xgr3=A1r:lx/4:A2r;%linspace(A1,A2,Nx);
Ygr3=B1r:ly/4:B2r;%linspace(B1,B2,Ny);
Xgr=sort([Xgr  Xgr3]);
Ygr=sort([Ygr  Ygr3]);
Xgr(find(abs(diff(Xgr))<1e-4)+1)=[];
Ygr(find(abs(diff(Ygr))<1e-4)+1)=[];
[XXgr YYgr]=meshgrid(Xgr,Ygr);

%forme identation

[~,aa]=min(abs(unique(XXgr)-b));
[~,bb]=min(abs(unique(XXgr)+b));
XXgr(:,aa)=b;
XXgr(:,bb)=-b;

[~,aa]=min(abs(unique(XXgr)-0));
[~,bb]=min(abs(unique(YYgr)-a));
YYgr(bb,aa)=a;

meshXY=[XXgr(:) YYgr(:)];
indxsup=find(XXgr(:)==0.);
indysup=find(YYgr(indxsup)==a);
% meshXY=[meshXY(1:indxsup(indysup)-1,:);0. 0.45;meshXY(indxsup(indysup):end,:)];
% indxsup=find(abs(meshXY(:,1)+lx)<1e-7);
% indysup=find(abs(meshXY(indxsup,2)-a)<1e-7);
% meshXY=[meshXY(1:indxsup(indysup),:);meshXY(indxsup(1),1) a+ly/2;meshXY(indxsup(indysup)+1:end,:)];
% indxsup=find(abs(meshXY(:,1)-lx)<1e-7);
% indysup=find(abs(meshXY(indxsup,2)-a)<1e-7);
% meshXY=[meshXY(1:indxsup(indysup),:);meshXY(indxsup(1),1) a+ly/2;meshXY(indxsup(indysup)+1:end,:)];
if strcmp(couplage,'join')
shapeindent=[A1 B2;-b 1;-b a;b a;b 1;A2 B2;A2 B1;A1 B1];
edgeindent=[(1:length(shapeindent)-1)' (2:length(shapeindent))'];
elseif strcmp(couplage,'zoom')
shapeindent=[A1 B2;-b B2;0 a; b B2;A2 B2;A2 B1;A1 B1];
edgeindent=[1 2;2 3;3 4;4 5;5 6;6 7;7 1];
end

%selection points de la grille dans la frome
[inindent,onindent]=inpoly(meshXY,shapeindent,edgeindent);


%maillage avec connec en plus
TRI2=DelaunayTri(meshXY([find(inindent);find(onindent)],1),meshXY([find(inindent);find(onindent)],2));

%derafinnement

[row1] = find(TRI2.X(:,2)<B1m - ly & TRI2.X(:,1)<A2m + lx & TRI2.X(:,1)>A1m -lx);
[row2] = find(TRI2.X(:,1)>A2m+lx & TRI2.X(:,2)<B2m+ly & TRI2.X(:,2)>B1m-ly);
[row3] = find(TRI2.X(:,1)<A1m-lx & TRI2.X(:,2)<B2m+ly & TRI2.X(:,2)>B1m-ly);
rowrem = [row1;row2;row3];
[rowcons1] = find((mod(TRI2.X(:,2),ly))==0);
%rowcons=[rowcons1;rowcons2];
rowrem2=setdiff(rowrem,rowcons1);
TRI2.X(rowrem2,:)=[];
[row1] = find(TRI2.X(:,2)<0.5-2*ly & TRI2.X(:,1)<2*lx & TRI2.X(:,1)>-2*lx);
[row2] = find(TRI2.X(:,1)>2*lx & TRI2.X(:,2)<a+2*ly & TRI2.X(:,2)>a-2*ly);
[row3] = find(TRI2.X(:,1)<-2*lx & TRI2.X(:,2)<a+2*ly & TRI2.X(:,2)>a-2*ly);
rowrem = [row1;row2;row3];
%valx = unique(TRI2.X(:,1));
%valrem = valx((abs((valx/lx)-round(valx/lx))>1e-6));
rowrem2 = find((abs((TRI2.X(:,1)/lx)-round(TRI2.X(:,1)/lx))>1e-6));
rowrem3 = intersect(rowrem,rowrem2);
TRI2.X(rowrem3,:)=[];
TRI=DelaunayTri(TRI2.X);

ptipass1=find(inpoly(TRI.X,[A1m-lx/4 B1m-lx/4;0 a]));
ptipass2=find(inpoly(TRI.X,[0 a;A2m+lx/4 B1m-lx/4]));
%ptipass = unique([ptipass1;ptipass2]);
lignpass1 = [ptipass1(1:end-1) ptipass1(2:end)];
lignpass2 = [ptipass2(1:end-1) ptipass2(2:end)];
TRI=DelaunayTri(TRI.X,[lignpass1;lignpass2]);

%ajustement forme indentation
meshXYindent=TRI.X;
connecXYindent=TRI.Triangulation;
ind1=find(abs(meshXYindent(:,1)+b)<1e-7);
ind2=find(abs(meshXYindent(:,1)-b)<1e-7);
rowxindent=min(ind1):max(ind2);
xindent=unique(meshXYindent(rowxindent,1));
for ijk=1:length(xindent)
    indtemp=find(abs(meshXYindent(:,1)-xindent(ijk))<1e-7);
    if strcmp(couplage,'zoom')
    meshXYindent(indtemp(end),2)=max(a-(1-a)/b*meshXYindent(indtemp(end),1),a+(1-a)/b*meshXYindent(indtemp(end),1));
    end
    contrain(ijk)=indtemp(end);
  %  [rowijk colindentconnec]=find(connecXYindent-indtemp(end)==0);
  %  connecXYindent(rowijk,:)=[];
end

connecmodif=connecXYindent;

if strcmp(couplage,'join')
[inindent,onindent]=inpoly(meshXYindent,[-b B2;-b a;b a;b B2],[1 2;2 3;3 4;4 1]);
contrain=find(onindent);
end
for ijk=1:length(contrain)    
    indzero=find(connecXYindent-contrain(ijk)==0);
    connecmodif(indzero)=0;
end
indremove=find(sum(connecmodif,2)==0);
connecXYindent(indremove,:)=[];
% for ijk=1:length(contrain)
%     for ijkl:length(contrain)
%         
%     end
% end

%TRindent = DelaunayTri(meshXYindent);

%connecXYindent=TRI.Triangulation;
%connecXYindent(connecremove)=[];

%[c rowindent1]=sort(abs(meshXYindent(:,2)-(a-(1-a)/b*meshXYindent(:,1))))
%[c rowindent2]=sort(abs(meshXYindent(:,2)-(a+(1-a)/b*meshXYindent(:,1))))


%TRIindent=DelaunayTri();

TRindent = TriRep(connecXYindent,meshXYindent );

%maillage micro

%Xgrm=A1m:lxfin:A2m;%linspace(A1,A2,Nx);
%Ygrm=B1m:lyfin:B2m;%linspace(B1,B2,Ny);
%[XXgrm YYgrm]=meshgrid(Xgrm,Ygrm);


%[~,aam]=min(abs(unique(XXgrm)-bm));
%[~,bbm]=min(abs(unique(XXgrm)+bm));
%XXgrm(:,aam)=bm;
%XXgrm(:,bbm)=-bm;

%[~,aam]=min(abs(unique(XXgrm)-0));
%[c,bbm]=min(abs(unique(YYgrm)-am));
%YYgrm(bbm,aam)=am;

%meshXYm=[XXgrm(:) YYgrm(:)];
shapeindentm=[A1m B2m;-bm B2m;0 am; bm B2m;A2m B2m;A2m B1m;A1m B1m];
edgeindent=[1 2;2 3;3 4;4 5;5 6;6 7;7 1];

%selection points de la grille dans la forme et contraintes
%[inindentm,onindentm]=inpoly(meshXYm,shapeindentm,edgeindent);
[inindent2]=inpoly(TRindent.X,shapeindentm,edgeindent);
%ETRindent = edges(TRindent);

meshXYemb = TRindent.X([find(inindent2)],:);
%meshXYemb = unique(meshXYemb,'rows','first');
%meshXYemb(find(abs(diff(meshXYemb(:,2)))<1e-7)+1,:)=[];
% 
% indemb1 = ismember(ETRindent(:,1),[find(inindent2);find(onindent2)]);
% indemb2 = ismember(ETRindent(:,2),[find(inindent2);find(onindent2)]);
% 
% ETRemb = ETRindent(indemb1 & indemb2,:);
% 

indemb1 = ismember(connecXYindent(:,1),find(inindent2));
indemb2 = ismember(connecXYindent(:,2),find(inindent2));
indemb3 = ismember(connecXYindent(:,3),find(inindent2));
connecXYemb=[];
for ijk=(find(indemb1 & indemb2 & indemb3))'
    V1 = TRindent.X(connecXYindent(ijk,1),:);
    V2 = TRindent.X(connecXYindent(ijk,2),:);
    V3 = TRindent.X(connecXYindent(ijk,3),:);
    indcons1 = find(ismember(meshXYemb,V1,'rows'));
    indcons2 = find(ismember(meshXYemb,V2,'rows'));
    indcons3 = find(ismember(meshXYemb,V3,'rows'));
    connecXYemb(end+1,:) = [indcons1 indcons2 indcons3];
end

%maillage avec connec en plus
TRItemp=TriRep(connecXYemb,meshXYemb);

[Xm Tm]=RefineMesh( TRItemp.X, TRItemp.Triangulation, nm );
TRindentm = TriRep(Tm,Xm);
% %ajustement forme indentation
% meshXYindentm=TRIm.X;
% connecXYindentm=TRIm.Triangulation;
% ind1m=find(meshXYindentm(:,1)==-bm);
% ind2m=find(meshXYindentm(:,1)==bm);
% rowxindentm=min(ind1m):max(ind2m);
% xindentm=unique(meshXYindentm(rowxindentm,1));
% for ijk=1:length(xindentm)
%     indtempm=find(meshXYindentm(:,1)==xindentm(ijk));
%     meshXYindentm(indtempm(end),2)=max(ap-(1-ap)/bp*meshXYindentm(indtempm(end),1),ap+(1-ap)/bp*meshXYindentm(indtempm(end),1));
%     contrainm(ijk)=indtempm(end);
%   %  [rowijk colindentconnec]=find(connecXYindent-indtemp(end)==0);
%   %  connecXYindent(rowijk,:)=[];
% end
% 
% connecmodifm=connecXYindentm;
% 
% for ijk=1:length(contrainm)    
%     [rawzerom,colzerom]=find(connecXYindentm-contrainm(ijk)==0);
%     for ijkl=1:length(rawzerom)
%     connecmodifm(rawzerom(ijkl),colzerom(ijkl))=0;
%     end
% end
% indremovem=find(sum(connecmodifm,2)==0);
% connecXYindentm(indremovem,:)=[];
% for ijk=1:length(contrain)
%     for ijkl:length(contrain)
%         
%     end
% end

%TRindent = DelaunayTri(meshXYindent);

%connecXYindent=TRI.Triangulation;
%connecXYindent(connecremove)=[];

%[c rowindent1]=sort(abs(meshXYindent(:,2)-(a-(1-a)/b*meshXYindent(:,1))))
%[c rowindent2]=sort(abs(meshXYindent(:,2)-(a+(1-a)/b*meshXYindent(:,1))))


%TRIindent=DelaunayTri();

%TRindentm = TriRep(connecXYindentm,meshXYindentm )

%TRIindent.X=meshXYindent;
%TRIindent.Triangulation=connecXYindent;
%TRIindent=DelaunayTri(meshXYindent(:,1),meshXYindent(:,2));

figure;triplot(TRindent)
xlim([A1-0.5 A2+0.5])
ylim([B1-0.5 B2+0.5])
axis equal
figure;triplot(TRindent,'color',[1 1 1]*.6, 'linewidth', 1)
xlim([A1-0.5 A2+0.5])
ylim([B1-0.5 B2+0.5])
axis equal
hold on
triplot(TRindentm,'color','k', 'linewidth', 0.1)
%plot(TRindent.X(:,1),TRindent.X(:,2),'r+')


if strcmp(couplage,'join')
LevelSet1=[-3 2;-bm 2;-bm am;bm am;bm 2;3 2;3 -1;-3 -1];
LevelSet2=[A1m B1m;A2m B1m;A2m 2;A1m 2];
elseif strcmp(couplage,'zoom')
    LevelSet1=[A1m B1m;A2m B1m;A2m B2m;A1m B2m];
    LevelSet2=[cx1 cy1;cx2 cy1;cx2 cy2;cx1 cy2];
   % LevelSet1=[-0.7 0.2;0.7 0.2;0.7 0.8;-0.7 0.8];%0.7 0.2 0.8
   % LevelSet2=[-0.5 0.4;0.5 0.4;0.5 0.6;-0.5 0.6];
end
hold on
plot(LevelSet1(:,1),LevelSet1(:,2),'r')
plot(LevelSet2(:,1),LevelSet2(:,2),'g')



%condition aux limites
shapeBC1=[A1 0;A1 1];
edgeBC1=[1 2];
[inindent1,onindent1]=inpoly(meshXYindent,shapeBC1,edgeBC1);

shapeBC2=[A2 0;A2 1];
edgeBC2=[1 2];
[inindent2,onindent2]=inpoly(meshXYindent,shapeBC2,edgeBC2);

strBC=repmat('U',1,length(find(onindent1))+length(find(onindent2)));
nodesBC=[find(onindent1)' find(onindent2)'];
valueBC=[-0.1*ones(1,length(find(onindent1))) 0.1*ones(1,length(find(onindent2)))];


%attribution parametre
Xmesomicro=sort(unique([-1:L/6:1]));
Ymesomicro=sort(unique([0:L/6:1]));
ind1=find(Xmesomicro==Xgr(1));
ind2=find(Xmesomicro==Xgr(end));
ind3=find(Ymesomicro==Ygr(1));
ind4=find(Ymesomicro==Ygr(end));
[kmesob kmicrob] = randomField3( 'lognormal', 'sinc2', L, 1, 0.1, Nmcmeso,Nmcmicro, Xmesomicro-A1, Ymesomicro-B1,[],dL);
[XXmesomicro YYmesomicro]=meshgrid(Xmesomicro,Ymesomicro);

% figure;pcolor(XXmesomicro(ind3:ind4,ind1:ind2)',YYmesomicro(ind3:ind4,ind1:ind2)',log(kmesob(ind1:ind2,ind3:ind4,end)));
% shading faceted
% colorbar
% figure;pcolor(XXmesomicro(ind3:ind4,ind1:ind2)',YYmesomicro(ind3:ind4,ind1:ind2)',log(kmicrob(ind1:ind2,ind3:ind4,end)));
% shading faceted
% colorbar
XXmesomicro=XXmesomicro';
YYmesomicro=YYmesomicro';

for ijk=1:Nmcmeso
    kter=kmesob(:,:,ijk);
   % Fmeso = TriScatteredInterp([reshape(XXmesomicro',length(Xmesomicro)*length(Ymesomicro),1) reshape(YYmesomicro',length(Xmesomicro)*length(Ymesomicro),1)],reshape(kter,length(Xmesomicro)*length(Ymesomicro),1));
   Fmeso = TriScatteredInterp(XXmesomicro(:),YYmesomicro(:),kter(:));
   Ftemp1=Fmeso(TRindent.X(:,1),TRindent.X(:,2));
    kmeso(:,:,ijk)=Ftemp1(TRindent.Triangulation);  %verif et interpolation 'TriScatteredInterp'
%kmeso(:,:,ijk)=ones(size(kmeso(:,:,ijk)));
end
%k = randomField2( 'lognormal', 'sinc2', 0.01, 1, 1, Nmc, Xgrm-A1m, Ygrm-B1m);
for ijk=1:Nmcmicro*Nmcmeso
    kter=kmicrob(:,:,ijk);
  %  Fmicro = TriScatteredInterp([reshape(XXmesomicro',length(Xmesomicro)*length(Ymesomicro),1) reshape(YYmesomicro',length(Xmesomicro)*length(Ymesomicro),1)],reshape(kter,length(Xmesomicro)*length(Ymesomicro),1));
   Fmicro = TriScatteredInterp(XXmesomicro(:),YYmesomicro(:),kter(:));
  Ftemp2=Fmicro(TRindentm.X(:,1),TRindentm.X(:,2));
    kmicro(:,:,ijk)=Ftemp2(TRindentm.Triangulation);
%kmicro(:,:,ijk)=ones(size(kmicro(:,:,ijk)));
end

Xmeso = TRindent.X(:,1);
Xmeso = Xmeso(TRindent.Triangulation);
Ymeso = TRindent.X(:,2);
Ymeso = Ymeso(TRindent.Triangulation);
Xmicro = TRindentm.X(:,1);
Xmicro = Xmicro(TRindentm.Triangulation);
Ymicro = TRindentm.X(:,2);
Ymicro = Ymicro(TRindentm.Triangulation);
figure;
patch(Xmeso',Ymeso',(kmeso(:,:,end)'));
hold on
shading faceted
patch(Xmicro',Ymicro',(kmicro(:,:,end)'));
shading faceted
xlim([A1-0.5 A2+0.5])
ylim([B1-0.5 B2+0.5])
colorbar

%    yselec = 0.5;
xselec = 0;
[xyselec pass] = sort(TRindent.X(abs(TRindent.X(:,1)-xselec)<1e-6,2));
    kyselec = Ftemp1((abs(TRindent.X(:,1)-xselec)<1e-6));
    kyselec = kyselec(pass);
    figure;plot(xyselec,kyselec);
    
    [xyselec pass] = sort(TRindentm.X(abs(TRindentm.X(:,1)-xselec)<1e-6,2));
    kyselec = Ftemp2((abs(TRindentm.X(:,1)-xselec)<1e-6));
    kyselec = kyselec(pass);
    hold on;plot(xyselec,kyselec);

% figure;
% patch(Xmeso',Ymeso',log(mean(mean(mean(kmeso,3)'))));
% hold on
% patch(Xmicro',Ymicro',log(mean(mean(mean(kmicro,3)'))));
% xlim([A1-0.5 A2+0.5])
% ylim([B1-0.5 B2+0.5])

support=1;

%LevelSet1=[];
%LevelSet2=[];
% cd Tests/
% if strcmp(couplage,'join')
% load('join2dstosto.mat');
% elseif strcmp(couplage,'zoom')
% load('zoom2dstosto.mat');
% end

solver='dmc';



model{1}.HomeFE.mesh.T=TRindent.Triangulation;
model{1}.HomeFE.mesh.X=TRindent.X;
model{1}.HomeFE.property=kmeso;
model{1}.HomeFE.load=zeros(size(model{1}.HomeFE.mesh.T));
model{1}.HomeFE.BC.type=strBC;
model{1}.HomeFE.BC.nodes=nodesBC;
model{1}.HomeFE.BC.value=valueBC;
model{1}.random.law='lognormal';
%model{1}.random.sig2=mean(mean(std(log(kmeso),[],1)));
%model{1}.random.min=0.9;
%model{1}.random.max=1.1;
model{1}.code='MonteCarloHomeFE';

model{2}.HomeFE.mesh.T=TRindentm.Triangulation;
model{2}.HomeFE.mesh.X=TRindentm.X;
model{2}.HomeFE.property=kmicro;
model{2}.HomeFE.load=zeros(size(model{2}.HomeFE.mesh.T));
model{2}.random.law='lognormal';
model{2}.HomeFE.BC=[];
%model{2}.random.sig2=mean(mean(std(log(kmicro),[],1)));
%model{2}.random.min=0.9;
%model{2}.random.max=1.1;
model{2}.code='MonteCarloHomeFE';

coupling{1}.levelSet = intersection( ...
            levelSet( 'square', LevelSet1(1,:), .25), ...
            levelSet( 'square',LevelSet2(1,:), .1, false) );
coupling{1}.cplval = {[1e-3 1-1e-3] [1e-3 1-1e-3]};
coupling{1}.freeval = {1e-3 1-1e-3};
coupling{1}.models=[1 2];
coupling{1}.mediator.support=support;
%coupling{1}.type=couplage;
coupling{1}.epsilon = 1e-3
coupling{1}.operator = 'H1';
coupling{1}.mediator.type='mesomicro';
%coupling{1}.LevelSet1.X=LevelSet1;
%coupling{1}.LevelSet1.T=  [[(1:length(LevelSet1)-1)' (2:length(LevelSet1))'];length(LevelSet1) 1];
%coupling{1}.LevelSet2.X=LevelSet2;
%coupling{1}.LevelSet2.T= [[(1:length(LevelSet2)-1)' (2:length(LevelSet2))'];length(LevelSet2) 1];
%cd ..
if strcmp(couplage,'join')
save('join2dstosto.mat','solver','coupling','model');
elseif strcmp(couplage,'zoom')
    save('zoom2dstosto.mat','solver','coupling','model');
end
%kmaj=k(sort([find(inindent) find(onindent)]));
% figure;plot(meshXYindent(:,1),meshXYindent(:,2),'+')
% xlim([-1.5 1.5])
% ylim([-0.5 1.5])

%maillage fin