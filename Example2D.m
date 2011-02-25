close all
clear all
% profile on

% load models
% load Tests/Coarse2D.mat
% load Tests/Medium2D.mat
% load Tests/Large2D.mat
load Tests/MonteCarloMedium2D.mat

% coupling{1}.mediator.support = 2;
model{2}.property = model{2}.property(:,:,1:50);

% start computation
[ sol, out ] = CArl( model, coupling, solver );

% profile viewer

% POSTTREATMENT
switch lower(solver)
    
% deterministic post treatment 2D
    case 'direct'
        X1 = model{1}.mesh.X; T1 = model{1}.mesh.T;
        X2 = model{2}.mesh.X; T2 = model{2}.mesh.T;
        figure; trimesh( T1, X1(:,1), X1(:,2), sol{1} );
        hold on ; trisurf( T2, X2(:,1), X2(:,2), sol{2}); shading flat

% stochastic post-treatment 2D
    case 'montecarlo'
        
        % mean displacements in 3D map
        X1 = model{1}.mesh.X; T1 = model{1}.mesh.T;
        X2 = model{2}.mesh.X; T2 = model{2}.mesh.T;
        figure; trimesh( T1, X1(:,1), X1(:,2), sol{1} );
        hold on ; trisurf( T2, X2(:,1), X2(:,2), sol{2}); shading flat

        % indices for line y=0
        ind1 = find( X1(:,2)==0 ); [x1,i1] = sort(X1(ind1,1)); ind1=ind1(i1);
        ind2 = find( X2(:,2)==0 ); [x2,i2] = sort(X2(ind2,1)); ind2=ind2(i2);

        % displacements
        u1 = sol{1}(ind1);
        u2 = out.u.MC.u(ind2,:);
        mu2 = mean(u2,2);
        su2 = std(u2,[],2);
        pmu2 = [ mu2+su2/sqrt(1-0.9); mu2(end:-1:1)-su2(end:-1:1)/sqrt(1-0.9) ];
        figure; fill( [x2;x2(end:-1:1)], pmu2, 'y' )
        hold on; plot( x1, u1, 'k-', x2, mu2, 'r--', x2, u2(:,1), 'b-' );

        s1 = diff( sol{1}(ind1) ) / mean(diff(x1)); s1 = [s1 s1]'; s1 = s1(:);
        xx1 = [x1(1:(end-1)) x1(2:end)]'; xx1 = xx1(:);
        xx2 = [x2(1:(end-1)) x2(2:end)]'; xx2 = xx2(:);
        s2 = zeros( length(xx2), size(out.u.MC.u(ind2,:),2));
        s2(1:2:end,:) = diff( out.u.MC.u(ind2,:), 1, 1 ) / mean(diff(x2));
        s2(2:2:end,:) = s2(1:2:end,:);
        ms2 = mean(s2,2);
        ss2 = std(s2,[],2);
        pms2 = [ ms2+ss2/sqrt(1-0.9); ms2(end:-1:1)-ss2(end:-1:1)/sqrt(1-0.9) ];
        figure; fill( [xx2;xx2(end:-1:1)], pms2, 'y' )
        hold on; plot( xx1, s1, 'k-', xx2, ms2, 'r--', xx2, s2(:,1), 'b-' );

end
