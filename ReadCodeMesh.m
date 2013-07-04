function model = ReadCodeMesh( model )
% READCODEMESH reads a mesh in a code-dependant format, and transforms it
% into a CArl format.
%
% syntax: mesh = ReadCodeMesh( model )
%
% accepted values for model.code are: 'HomeFe', 'MonteCarloHomeFE', 'FE2D',
% 'Comsol' and 'Beam'

% creation: R. Cottereau 01/2010

% switch on the type of code used
switch model.code
    
    % HOMEFE
    case {'HomeFE','MonteCarloHomeFE','FE2D'}
        model.mesh = TRI6( model.HomeFE.mesh.T, model.HomeFE.mesh.X, false );

    % COMSOL
    case 'Comsol'
        wd = pwd;
        cd( model.meshpath );
        eval( ['mesh = ' model.meshfile ';' ]);
        X = mesh.mesh('mesh1').getVertex';
        T = double( mesh.mesh('mesh1').getElem('tri')' +1 );
        model.mesh = TRI6( T, X, false );
        eval( ['cd ' wd ';' ]);
    
    % TIMOSCHENKO BEAM CODE
    case  'Beam'
        model.mesh = model.Beam.mesh;
        L = model.Beam.L;
        Xb = model.Beam.mesh.X;
        if size(Xb,2)>size(Xb,1)
            Xb = Xb';
        end
        if isfield( model, 'mesh' )
            Y = [-L/2 -L/4 0 L/4 L/2];
        else
            Y = [-3*L/4 -L/2 -L/4 0 L/4 L/2 3*L/4];
        end
        [X,Y] = meshgrid( Xb, Y );
        XY = [ X(:) Y(:) ];
        dt = DelaunayTri( XY );
        model.virtualMesh2D = TRI6( dt.Triangulation, dt.X, false );
        [~,model.v2mX] = ismember( dt.X(:,1), Xb );
        model.v2mT = zeros(model.virtualMesh2D.Ne,1);
        XindT = model.v2mX(dt.Triangulation);
        for i1 = 1:size(model.mesh.T,1)
            model.v2mT( all( XindT==model.mesh.T(i1,1) | ...
                             XindT==model.mesh.T(i1,2), 2 ) ) = i1;
        end
        keyboard
    % error
    otherwise
        error( 'unknown code' );
end
        