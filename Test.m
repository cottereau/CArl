function Test(type)
% launch a test or a series of tests
% The existing tests are the following
%
% ONE-DIMENSIONAL TESTS
% 'zoom1D' for a zoom case in 1D (deterministic continuum-continuum)
% 'join1D' for a junction case in 1D (deterministic continuum-continuum)
% 'force1D'
% 'MC1D' for a joint case in 1D (stochastic continuum-continuum)
% 'MC1Dstosto' for a joint case in 1D (stochastic-stochastic continuum)
% 'MC1D_BC_u10'
% 'MC1Dstosto2' for a joint case in 1D (stochastic-stochastic coupling)
% 'zoom1Dstosto' for a zoom case in 1D (stochastic-stochastic coupling)
%
% TWO-DIMENSIONAL TESTS
% 'join2D' for a junction case in 2D (deterministic continuum-continuum)
% 'zoom2D'
% 'force2D'
% 'comsol2D'
% 'nonembedded2D_1'
% 'zoom2Dindent'
%
% the possible series are
% 'short' for only the short tests
% 'comsol' for the comsol tests

% R. Cottereau 05/2010
% TO DO
% - a test with three models and two couplings
%
% modif YLG 02/2013: add displacement plots for det-sto coupling
% modif YLG 03/2013: add the meso-micro coupling

% default
if nargin==0
    type = 'short';
end

% selection of type of test
switch lower(type)
    
    case {'zoom1d', 'force1d', 'join1d'}
        load(['Tests/' type '.mat']);
        sol = CArl( model, coupling, solver );
        plottest1D( model, sol, ref );
        
    case {'join2d', 'force2d', 'zoom2d', 'zoom2dstosto', 'comsol2d', ...
          'nonembedded2d_1', 'zoom2dindent'}
        load(['Tests/' type '.mat']);
        [ sol, out ] = CArl( model, coupling, solver );
        plottest2D( out.model, sol,out );

    case {'mc1d','mc1d_bc_u10','mc1dstosto2','zoom1dstosto'}
        load(['Tests/' type '.mat']);
        [ sol, out ] = CArl( model, coupling, solver );
        plotteststochastic2( model, sol, out, ref );
        
    case 'short'
        Test('zoom1D');
        Test('join1D');
        Test('force1D');
        Test('MC1D');
        Test('MC1D_BC_u10') ;
        Test('mc1dstosto2');
        Test('zoom1Dstosto');
        Test('2d');
        
    case 'comsol'
        Test('comsol2D');
        
    case '2d'
        Test('join2D');
        Test('zoom2D');
        Test('zoom2Dindent');
        Test('NonEmbedded2D_1');
        Test('force2D');
        
    otherwise
        error('unknown test case')

end

% FUNCTION PLOTTEST 1D
function plottest1D( model, sol, ref )
x1 = model{1}.mesh.X( model{1}.mesh.T )';
x2 = model{2}.mesh.X( model{2}.mesh.T )';
u1 = sol{1}';
u2 = sol{2}';
figure; plot( x1, u1, 'bx-', x2, u2, 'ro--' )
if ~isempty(ref)
    hold on; plot(ref.X,ref.u,'k--')
end
s1 = diff(u1) ./ diff(x1); s1 = [s1; s1];
s2 = diff(u2) ./ diff(x2); s2 = [s2; s2];
figure; plot( x1, s1, 'bx-', x2(1,:), s2(1,:), 'r' );
if ~isempty(ref)
    xref = zeros( 2*(length(ref.X)-1), 1);
    xref(1:2:end,:) = ref.X(1:end-1) ;
    xref(2:2:end,:) = ref.X(2:end) ;
    sref = zeros( 2*(length(ref.X)-1), 1);
    sref(1:2:end,:) = diff( ref.u, 1 ) / mean(diff(ref.X));
    sref(2:2:end,:) = sref(1:2:end,:);
    hold on; plot(xref,sref,'k--')
end

% FUNCTION PLOTTEST 2D
function plottest2D( model, sol,out )
T1 = model{1}.mesh.T3;
T2 = model{2}.mesh.T3;
X1 = model{1}.mesh.X3;
X2 = model{2}.mesh.X3;
figure; trimesh( T1, X1(:,1), X1(:,2), sol{1} ); 
hold on; trisurf( T2, X2(:,1), X2(:,2), sol{2} ); 
colorbar; view(3)

        if strcmp(model{1}.code,'MonteCarloHomeFE')

T1 = out.model{1}.mesh.tri3.Triangulation;
        X1 = out.model{1}.mesh.tri3.X;
        x1 = X1(:,1); y1 = X1(:,2);
        Xc1 = [mean(x1(T1),2) mean(y1(T1),2)];
        uMC1 = out.model{1}.uMC;
        alpha1 = out.model{1}.alpha;
        T2 = out.model{2}.mesh.tri3.Triangulation;
        X2 = out.model{2}.mesh.tri3.X;
        x2 = X2(:,1); y2 = X2(:,2);
        Xc2 = [mean(x2(T2),2) mean(y2(T2),2)];
        uMC2 = out.model{2}.uMC;
        alpha2 = out.model{2}.alpha;
        Nmc = size(model{2}.HomeFE.property,3);
        Nmcmeso = size(model{1}.HomeFE.property,3);
        
        % plot mesh
        figure; trimesh(T1,X1(:,1),X1(:,2),'color',[1 1 1]*.6, 'linewidth', 1)
        hold on; trimesh(T2,X2(:,1),X2(:,2),'color','k', 'linewidth', 0.1)
        set( gca, 'PlotBoxAspectRatio', [2 1 1] ); box off;
        xlabel( 'Position X [m] '); ylabel( 'Position Y [m] ');
        cplx = [-.125 -.125 .125 .125 .05 .05 -.05 -.05 -.125];
        cply = [ .625 .375 .375 .625 .55 .45 .45 .55 .625];
        figure; trimesh(T1,X1(:,1),X1(:,2),'color',[1 1 1]*.6, 'linewidth', 1)
        hold on; trimesh(T2,X2(:,1),X2(:,2),'color','k', 'linewidth', 0.5)
        hold on; plot( cplx, cply, 'k-', 'linewidth', 3 )
        set( gca, 'PlotBoxAspectRatio', [.42 .41 1] ); box off;
        xlabel( 'Position X [m] '); ylabel( 'Position Y [m] ');
        set( gca, 'XLim', [-1 1]*.21, 'YLim', [.29 .7] )
        
        
        sol1 = sol{1};
        sol2 = sol{2};
        figure; patch(x1(T1)',y1(T1)',sol1(T1)',sol1(T1)')
hold on; patch(x2(T2)',y2(T2)',sol2(T2)',sol2(T2)')
view(3)
        % plot gr(alpham*E(wm)+alphaM*wM) on the overall domain
                
        for ijk = 1:Nmcmeso
            wmtemp = mean(uMC2(:,(ijk-1)*Nmc/Nmcmeso+1:ijk*Nmc/Nmcmeso),2);
            grwmtemp = trigradient(T2,x2,y2,wmtemp);
            wmtemp = wmtemp(T2);
            wm(:,:,ijk) = full(wmtemp);
            grwm(:,:,ijk)=grwmtemp(T2);
            FwM = TriScatteredInterp(x1,y1,full(uMC1(:,ijk*Nmc/Nmcmeso)));
            wMX2temp = FwM(X2);
            wMX2(:,:,ijk) = wMX2temp(T2);
            grwMX2temp = trigradient(T2,x2,y2,wMX2temp);
            grwMX2(:,:,ijk)=grwMX2temp(T2);
        end
        alphamX2 = interp(alpha2,X2);
        gralphamX2 = trigradient(T2,x2,y2,alphamX2);
        amwm = repmat(alphamX2(T2),[1 1 Nmcmeso]).*wm;
        alphaMX2 = interp(alpha1,X2);
        gralphaMX2 = trigradient(T2,x2,y2,alphaMX2);
        
        graMwMamEamX2 = repmat(alphamX2(T2),[1 1 Nmcmeso]).*grwm + repmat(alphaMX2(T2),[1 1 Nmcmeso]).*grwMX2 + repmat(gralphamX2(T2),[1 1 Nmcmeso]).*wm + repmat(gralphaMX2(T2),[1 1 Nmcmeso]).*wMX2;
        aMwMamEamX2 = repmat(alphamX2(T2),[1 1 Nmcmeso]).*wm + repmat(alphaMX2(T2),[1 1 Nmcmeso]).*wMX2;
        
        Xbf2=freeBoundary(out.model{2}.mesh);
        [outdomainmeso ondomainmeso]=(inpoly(X1,Xbf2.X{1},Xbf2.T{1}));
        for ijk=1:Nmcmeso
        aMwMmesotemp = uMC1(:,ijk*Nmc/Nmcmeso);
       % aMwMmeso(outdomainmeso)=NaN;
        graMwMmesotemp = trigradient(T1,x1,y1,aMwMmesotemp);
        graMwMmesotemp(outdomainmeso)=0;
        graMwMmeso(:,:,ijk)=graMwMmesotemp(T1);
        aMwMmeso(:,:,ijk) = full(aMwMmesotemp(T1));
        end
   for ijk = 1:Nmc
       wmtemp = uMC2(:,ijk);
       gramwmtemp = trigradient(T2,x2,y2,wmtemp);
        gramwm(:,:,ijk) = gralphamX2(T2).*wmtemp(T2) +alphamX2(T2).*gramwmtemp(T2);
   end
        figure;
        patch(x2(T2)',y2(T2)',mean(gramwm,3)');shading interp
                
        figure;
        patch(x2(T2)',y2(T2)',var(gramwm,[],3)');shading interp
        
                figure;
        patch(x1(T1)',y1(T1)',mean(aMwMmeso,3)');shading interp
                        hold on
       patch(x2(T2)',y2(T2)',mean(aMwMamEamX2,3)',mean(aMwMamEamX2,3)');shading interp
        colorbar
        
        
        figure;
        patch(x1(T1)',y1(T1)',mean(graMwMmeso,3)');shading interp
       patch(x2(T2)',y2(T2)',mean(graMwMamEamX2,3)');shading interp
        colorbar
                hold on
         
xselec = 0;
[rowxselec colxselec] = find([x1(T1);x2(T2)]==xselec);
yxselecb = [y1(T1);y2(T2)];
uMarl = cat(1,aMwMmeso,aMwMamEamX2);
gruMarl = cat(1,graMwMmeso,graMwMamEamX2);
for ijk=1:length(rowxselec)
yxselec(ijk) = yxselecb(rowxselec(ijk),colxselec(ijk));
uMarlselec(ijk,:) =uMarl(rowxselec(ijk),colxselec(ijk),:);
gruMarlselec(ijk,:)=gruMarl(rowxselec(ijk),colxselec(ijk),:);
end
[yxselecdef m n] = unique(yxselec,'last');
uMarlselecdef = uMarlselec(m,:); 
gruMarlselecdef = gruMarlselec(m,:);

pc=0.9;

yxselecfill = [yxselecdef(1:end-1);yxselecdef(2:end)];yxselecfill=yxselecfill(:);
yxselecfill = [yxselecfill;yxselecfill(end:-1:1)];
meanuMarl=transpose(mean(cat(2,reshape(uMarlselecdef(1:end-1,:),[length(yxselecdef)-1,1,Nmcmeso]),reshape(uMarlselecdef(2:end,:),[length(yxselecdef)-1,1,Nmcmeso])),3));
meanuMarl=meanuMarl(:);
stduMarl=transpose(std(cat(2,reshape(uMarlselecdef(1:end-1,:),[length(yxselecdef)-1,1,Nmcmeso]),reshape(uMarlselecdef(2:end,:),[length(yxselecdef)-1,1,Nmcmeso])),[],3));
stduMarl=stduMarl(:);
probuMarl=[meanuMarl+stduMarl/sqrt(1-pc);(meanuMarl(end:-1:1)-stduMarl(end:-1:1)/sqrt(1-pc))];

figure;
fill( yxselecfill, probuMarl, 'y' );
hold on;
plot(yxselecdef,uMarlselecdef(:,end),'r');
plot(yxselecdef,mean(uMarlselecdef,2));
title('uMarl x=0')
                
 
meangruMarl=transpose(mean(cat(2,reshape(gruMarlselecdef(1:end-1,:),[length(yxselecdef)-1,1,Nmcmeso]),reshape(gruMarlselecdef(2:end,:),[length(yxselecdef)-1,1,Nmcmeso])),3));
meangruMarl=(meangruMarl(:));
stdgruMarl=transpose(std(cat(2,reshape(gruMarlselecdef(1:end-1,:),[length(yxselecdef)-1,1,Nmcmeso]),reshape(gruMarlselecdef(2:end,:),[length(yxselecdef)-1,1,Nmcmeso])),[],3));
stdgruMarl=stdgruMarl(:);
probgruMarl=[meangruMarl+stdgruMarl/sqrt(1-pc);(meangruMarl(end:-1:1)-stdgruMarl(end:-1:1)/sqrt(1-pc))];

figure;
fill( yxselecfill, probgruMarl, 'y' );
hold on;
plot(yxselecdef,gruMarlselecdef(:,end),'r');
plot(yxselecdef,mean(gruMarlselecdef,2));
title('gruMarl x=0')      
end


% FUNCTION PLOTTESTSTOCHASTIC
function plotteststochastic2( model, ~, out, ref )
pc = 0.9;
x1 = model{1}.mesh.X( model{1}.mesh.T )';
x2 = model{2}.mesh.X( model{2}.mesh.T )';
dx2 = diff(x2);
dx1 = diff(x1);
u1=mean(out.model{1}.auMC,3)';uu1=u1(:);
stdu1=std(out.model{1}.auMC,0,3)';stdu1=stdu1(:);
probu1=[uu1+stdu1/sqrt(1-pc);(uu1(end:-1:1)-stdu1(end:-1:1)/sqrt(1-pc))];
xxx1 = [ (x1(:));(x1(end:-1:1))'];
figure;
if (strcmp(out.coupling{1}.mediator.type,'mesomicro'))
    fill( xxx1, probu1, 'y' )
end
hold on
%u1b = (sol{1})';
%u2b=(sol{2})';
u2=mean(out.model{2}.auMC,3)';uu2=u2(:);%mean u2
stdu2=std(out.model{2}.auMC,0,3)';stdu2=stdu2(:);%std u2
probu2=[uu2+stdu2/sqrt(1-pc);uu2(end:-1:1)-stdu2(end:-1:1)/sqrt(1-pc)];
xxx2 = [ x2(:);x2(end:-1:1)'];
fill( xxx2(:), probu2(:), 'y' )
plot(x1,u1,'bx-',x2,u2,'r--');


%s1b = diff(u1) ./ diff(x1);% s1b = [s1b; s1b];
%s2b = diff(u2) ./ diff(x2);% s2b = [s2b; s2b];

s1 = squeeze(diff(out.model{1}.auMC,[],2)) ...
           ./ repmat(dx1',[1 size(out.model{1}.auMC,3)]);
ms1 = mean(s1,2); ms1 = [ms1 ms1]'; ms1 = ms1(:);
ss1 = std(s1,[],2); ss1 = [ss1 ss1]'; ss1 = ss1(:);
pms1 = [ ms1+ss1/sqrt(1-pc); ms1(end:-1:1)-ss1(end:-1:1)/sqrt(1-pc) ];
xx1 = x1(:); xx1 = [ xx1; xx1(end:-1:1) ];
figure;
if (strcmp(out.coupling{1}.mediator.type,'mesomicro'))
fill( xx1, pms1, 'y' )
end
hold on;
s2 = squeeze(diff(out.model{2}.auMC,[],2)) ...
           ./ repmat(dx2',[1 size(out.model{2}.auMC,3)]);
ms2 = mean(s2,2); ms2 = [ms2 ms2]'; ms2 = ms2(:);
ss2 = std(s2,[],2); ss2 = [ss2 ss2]'; ss2 = ss2(:);
pms2 = [ ms2+ss2/sqrt(1-pc); ms2(end:-1:1)-ss2(end:-1:1)/sqrt(1-pc) ];
xx2 = x2(:); xx2 = [ xx2; xx2(end:-1:1) ];
fill( xx2, pms2, 'y' )
plot( x1(:), ms1, 'bx-', x2(:), ms2, 'ro--' );
if ~isempty(ref)
    xref = [ref.X(1:(end-1)) ref.X(2:end)]'; xref = xref(:);
    sref = diff(ref.u) ./ diff(ref.X); sref = [sref; sref];
    hold on; plot(xref,sref,'k--')
end

