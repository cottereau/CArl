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
    case {'homefe'}
        T = model.HomeFE.mesh.T;
        X = model.HomeFE.mesh.X;
        d = size(X,2);
        if d==1
            model.mesh = INT3( T, X, false );
        elseif d==2
            model.mesh = TRI6( T, X, false );
        else
            error('case not treated yet')
        end

    % FE2D elastic code
    case 'fe2d'
        model.mesh = TRI6( model.FE2D.T, model.FE2D.X, false );

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
    case  'beam'
        if isfield( model.Beam, 'dir' )
            dir = model.Beam.dir;
        else
            dir =[];
        end
        if isfield( model.Beam, 'x0' )
            x0 = model.Beam.x0;
        else
            x0 = [];
        end
        model.mesh = INT3( model.Beam.T, model.Beam.X, false, dir, x0 );
        
    % error
    otherwise
        error( 'unknown code' );
end
        