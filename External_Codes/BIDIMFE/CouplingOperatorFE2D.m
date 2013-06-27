function [ x, y, C ] = CouplingOperatorFE2D( operator, mesh, opt )


[Nne,ne] = size(mesh.Triangulation);

% computation of MassMatrix
m = struct( 'mesh', mesh, ...
    'property', ones(Nne,ne), ...
    'load', zeros(Nne,ne,2), ...
    'BC', [] );

elements = m.mesh.Triangulation;
coordinates = m.mesh.X;

[ Mass1,Mass2 ] = mass_coupling_Timo(elements,coordinates);
[x2,y2] = find(Mass2);
indtemp = (find(Mass2));
C = Mass2(indtemp);
[x,y] = find(Mass1);
maxx = size(Mass1,1);
x = [x;x2];
y = [y;y2+maxx];
indtemp = (find(Mass1));
C = [Mass1(indtemp);C];
% choice of the coupling operator
switch operator
    
    % L2 coupling
    case 'L2'
        
        % H1 coupling
    case 'H1'
        [stiffC1,stiffC2] = stifness_coupling_Timo( elements,coordinates);
        [z2,k2] = find(stiffC2);
        indtemp= find(stiffC2);
        K2=(stiffC2(indtemp));
        x = [ z2; x ];
        y = [ k2+maxx; y ];
        C = [ opt.kappa*K2; C ];
        
        [z,k] = find(stiffC1);
        indtemp= find(stiffC1);
        K=(stiffC1(indtemp));
        x = [ z; x ];
        y = [ k; y ];
        C = [ opt.kappa*K; C ];
        % unknown coupling operator
    otherwise
        error('unknown coupling operator')
end