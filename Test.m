function Test(type)
% launch series of tests
% two possibilites: Test('full') to launch all test, which may be long,
% or Test('short') (or simply Test) which will only launch the fast ones

% R. Cottereau 05/2010

% default
if nargin==0
    type = 'short';
end

% FAST TESTS
load Tests/Coarse1D.mat;
sol = CArl( model, coupling, solver );
plottest1D( model, sol );

load Tests/Medium1D.mat;
sol = CArl( model, coupling, solver );
plottest1D( model, sol );

load Tests/StochasticCoarse1D.mat;
model{2}.property = model{2}.property(:,1:10);
[ sol, out ] = CArl( model, coupling, solver );
plotteststochastic( model, sol, out );

load Tests/StochasticMedium1D.mat;
model{2}.property = model{2}.property(:,1:10);
[ sol, out ] = CArl( model, coupling, solver );
plotteststochastic( model, sol, out );

load Tests/Coarse2D.mat;
sol = CArl( model, coupling, solver );
plottest2D( model, sol );


% SLOW TESTS
switch type
    case 'full'
        load Tests/Fine1D.mat;
        sol = CArl( model, coupling, solver );
        plottest1D( model, sol );

        load Tests/StochasticCoarse1D.mat;
        [ sol, out ] = CArl( model, coupling, solver );
        plotteststochastic( model, sol, out );

        load Tests/StochasticMedium1D.mat;
        [ sol, out ] = CArl( model, coupling, solver );
        plotteststochastic( model, sol, out );


end


% FUNCTION PLOTTEST 1D
function plottest1D( model, sol )
figure; plot(model{1}.mesh.X,sol{1}, 'b', model{2}.mesh.X,sol{2},'r')

% FUNCTION PLOTTEST 2D
function plottest2D( model, sol )
figure;
trimesh( model{1}.mesh.T, model{1}.mesh.X(:,1), model{1}.mesh.X(:,2), sol{1});
hold on ;
trisurf( model{2}.mesh.T, model{2}.mesh.X(:,1), model{2}.mesh.X(:,2), sol{2});


% FUNCTION PLOTTESTSTOCHASTIC
function plotteststochastic( model, sol, out )
X1 = model{1}.mesh.X;
s1 = diff( sol{1} ) / mean(diff(X1)); s1 = [s1 s1]'; s1 = s1(:);
x1 = [X1(1:(end-1)) X1(2:end)]'; x1 = x1(:);
X2 = model{2}.mesh.X;
x2 = [X2(1:(end-1)) X2(2:end)]'; x2 = x2(:);
s2 = zeros( length(x2), size(out.u.MC.u,2));
s2(1:2:end,:) = diff( out.u.MC.u, 1, 1 ) / mean(diff(X2));
s2(2:2:end,:) = s2(1:2:end,:);
ms2 = mean(s2,2);
ss2 = std(s2,[],2);
pms2 = [ ms2+ss2/sqrt(1-0.9); ms2(end:-1:1)-ss2(end:-1:1)/sqrt(1-0.9) ];
figure; fill( [x2;x2(end:-1:1)], pms2, 'y' )
hold on; plot( x1, s1, 'k-', x2, ms2, 'r--' );

