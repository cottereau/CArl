function model = ReadCodeMesh( model )
% READCODEMESH reads a mesh in a code-dependant format, and transforms it
% into a CArl format.
%
% syntax: mesh = ReadCodeMesh( model )
%
%  model is a structured array that should contain a field 'code'
%        indicating the exterior code that is being used. Other fields
%        are code-dependent. Fields 'mesh' should not be used because they
%        used internally by CArl.
%
% accepted values for model.code are: 'HomeFE', 'MonteCarloHomeFE', 'FE2D',
% 'Comsol' and 'Beam'

% creation: R. Cottereau 01/2010

% switch on the type of code used
switch lower(model.code)
    
    % HOMEFE
    case {'homefe','montecarlohomefe','fe2d'}
        model.mesh = TRI6( model.HomeFE.mesh.T, model.HomeFE.mesh.X, false );

    % COMSOL
    case 'comsol'
        wd = pwd;
        cd( model.meshpath );
        eval( ['mesh = ' model.meshfile ';' ]);
        X = mesh.mesh('mesh1').getVertex';
        T = double( mesh.mesh('mesh1').getElem('tri')' +1 );
        model.mesh = TRI6( T, X, false );
        eval( ['cd ' wd ';' ]);
    
    % TIMOSCHENKO BEAM CODE
    % we are assuming here that the local referential of the beam is along 
    % for a line in direction [1 0] and going through point [0 0]
    case  'beam'
        % default values
        Xb = model.Beam.mesh.X;
        if size(Xb,2)>size(Xb,1)
            Xb = Xb';
        end
        if isfield( model.Beam, 'L' )
            L = model.Beam.L;
        else
            L = 1;
        end
        if isfield( model.Beam, 'X0' )
            X0 = model.Beam.X0;
        else
            X0 = [0 0];
        end
        if isfield( model.Beam, 'dir' )
            dir = model.Beam.dir;
        else
            dir = [1 0];
        end
        % create 1D mesh of INT3
        model.mesh = INT3( model.Beam.mesh.T, Xb );
        % create virtual 2D mesh of TRI6 and transformation basis
        orth = [-dir(2) dir(1)];
        Y = [-L/2 0 L/2];
        [X,Y] = meshgrid( Xb, Y );
        XY = X(:)*dir + Y(:)*orth;
        XY = [ XY(:,1)+X0(1)  XY(:,2)+X0(2) ];
        dt = DelaunayTri( XY );
        model.virtualMesh2D = TRI6( dt.Triangulation, dt.X, false );
        Nl = length(Xb);
        model.v2mX = reshape( repmat( 1:2:2*Nl, [3 1] ), 3*Nl, 1 );
        model.v2mT = zeros( model.virtualMesh2D.Ne, 1 );
        Tloc = model.v2mX(model.virtualMesh2D.T3);
        for i1 = 1:model.mesh.Ne
            ind = all( ismember(Tloc,model.mesh.T(i1,:)), 2 );
            model.v2mT(ind) = i1;
        end
    % error
    otherwise
        error( 'unknown code' );
end
        