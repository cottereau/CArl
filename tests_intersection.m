function tests_intersection(cas)

geps = 1e-10;
close all

if nargin<1
    cas = 'long';
end

switch cas
    case 'short'
        
        % mesh 1
        node1 = [0 0;1 0;1 1;0 1]; [X1,T1] = mesh2d(node1);
        node2 = node1 + .25; [X2,T2] = mesh2d(node2);
        [Ti,Xi] = intersect_trimesh( T1, X1, T2, X2, geps );
        close(1:2)
        plotmeshes( 1, T1, X1, T2, X2, Ti, Xi )

        % mesh 2
        [X2,T2] = mesh2d(node2,[],struct('hmax',0.1));
        [Ti,Xi] = intersect_trimesh( T1, X1, T2, X2, geps );
        close(2)
        plotmeshes( 2, T1, X1, T2, X2, Ti, Xi )

        % mesh 3
        node1 = [0 0;1 0;1 1;0 1]; [X1,T1] = mesh2d(node1);
        node2 = [.2 0;1.2 0.2;1 1.2;0 1]; [X2,T2] = mesh2d(node2);
        [Ti,Xi] = intersect_trimesh( T1, X1, T2, X2, geps );
        close(3:4)
        plotmeshes( 3, T1, X1, T2, X2, Ti, Xi )
        
        % mesh 4
        [X2,T2] = mesh2d(node2,[],struct('hmax',0.1));
        [Ti,Xi] = intersect_trimesh( T1, X1, T2, X2, geps );
        close(4)
        plotmeshes( 4, T1, X1, T2, X2, Ti, Xi )

    case 'long'
        
        % mesh 1
        node1 = [0 0;1 0;1 1;0 1];
        node2 = node1 + .25;
        [X1,T1] = mesh2d(node1,[],struct('hmax',0.1));
        [X2,T2] = mesh2d(node2,[],struct('hmax',0.1));
        [Ti,Xi] = intersect_trimesh( T1, X1, T2, X2, geps );
        close(1:2)
        plotmeshes( 1, T1, X1, T2, X2, Ti, Xi )

        % mesh 2
        node2 = [.2 0;1.2 0.2;1 1.2;0 1];
        [X2,T2] = mesh2d(node2,[],struct('hmax',0.1));
        [Ti,Xi] = intersect_trimesh( T1, X1, T2, X2, geps );
        close(2)
        plotmeshes( 2, T1, X1, T2, X2, Ti, Xi )

end

% =============================================
function plotmeshes( h, T1, X1, T2, X2, Ti, Xi )
figure(h); trimesh(T1,X1(:,1),X1(:,2),'color','k');
hold on; trimesh(T2,X2(:,1),X2(:,2),'color','r');
hold on; trimesh(Ti,Xi(:,1),Xi(:,2),'color','g');