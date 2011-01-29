close all
clear all
% profile on

% load models
% load Tests/Coarse1D.mat
% load Tests/Medium1D.mat
% load Tests/Large1D.mat
% load Tests/MonteCarloCoarse1D.mat
load Tests/MonteCarloMedium1D.mat

% coupling{1}.mediator.support = 2;
model{2}.property = model{2}.property(:,:,1:100);

% start computation
[ sol, out ] = CArl( model, coupling, solver );

% profile viewer

% POSTTREATMENT
switch solver
    
% deterministic post treatment 1D
case 'direct'
    figure; plot( model{1}.mesh.X, sol{1}, 'b' );
    hold on; plot( model{2}.mesh.X, sol{2}, 'r' );

% stochastic post-treatment 1D
case 'MonteCarlo'
    % positions
    X1 = model{1}.mesh.X; dX1 = mean(diff(X1));
    x1 = [X1(1:(end-1)) X1(2:end)]'; x1 = x1(:);
    X2 = model{2}.mesh.X; dX2 = mean(diff(X2)); NX2 = length(X2);
    x2 = [X2(1:(end-1)) X2(2:end)]'; x2 = x2(:); Nx2 = length(x2);

    % displacements
    u1 = sol{1};
    u2 = out.u.MC.u(1:NX2,:);
    mu2 = mean(u2,2);
    su2 = std(u2,[],2);
    pmu2 = [ mu2+su2/sqrt(1-0.9); mu2(end:-1:1)-su2(end:-1:1)/sqrt(1-0.9) ];
    figure; fill( [X2;X2(end:-1:1)], pmu2, 'y' )
    hold on; plot( X1, u1, 'k-', X2, mu2, 'r--', X2, u2(:,1), 'b-' );
    
    % gradients
    s1 = diff( sol{1} ) / dX1; s1 = [s1 s1]'; s1 = s1(:);
    s2 = zeros( Nx2, size(out.u.MC.u,2));
    s2(1:2:end,:) = diff( out.u.MC.u(1:NX2,:), 1, 1 ) / dX2;
    s2(2:2:end,:) = s2(1:2:end,:);
    ms2 = mean(s2,2);
    ss2 = std(s2,[],2);
    pms2 = [ ms2+ss2/sqrt(1-0.9); ms2(end:-1:1)-ss2(end:-1:1)/sqrt(1-0.9) ];
    figure; fill( [x2;x2(end:-1:1)], pms2, 'y' )
    hold on; plot( x1, s1, 'k-', x2, ms2, 'r--', x2, s2(:,1), 'b-' );

end

