function model = ReadCodeMesh( model )
% READCODEMESH reads a mesh in a code-dependant format, and transforms it
% into a CArl format.
%
% syntax: mesh = ReadCodeMesh( model )
%
% accepted values for model.code are: 'HomeFe', 'MonteCarloHomeFE', and
% 'Comsol'
%
% 1) 'HomeFe', 'MonteCarloHomeFE'
%    the field model.mesh.X and model.mesh.T should exist and be:
%      -X: coordinates matrix [Nn*d double matrix], with de the space
%          dimension
%      -T: connectivity matrix (only simplices) [Ne*(d+1) integer matrix]
%
% 2) 'Comsol'
%    the field model.meshpath and model.meshfile should exist and be:
%      -meshpath: path to meshfile
%      -meshfile: name of a matlab command to read the mesh information
%      NB: the mesh should be stored as 'mesh1'

% R. Cottereau 01/2010

% switch on the type of code used
switch model.code
    
    % HOMEFE
    case {'HomeFE','MonteCarloHomeFE','FE2D'}
        model.mesh = TRI6( model.HomeFE.mesh.T, model.HomeFE.mesh.X, false );
    case  {'TimobeamFE'}
        tmpX = model.HomeFE.mesh.X;
        tmpY = [-model.HomeFE.L/2 0 model.HomeFE.L/2];
        [tmpXX tmpYY] = meshgrid(tmpX,tmpY);
        coordinates = [tmpXX(:) tmpYY(:)];
        tmpmesh = DelaunayTri(coordinates);
        elements3 = tmpmesh.Triangulation;
        model.mesh = TRI6( elements3, coordinates, false );
    % COMSOL
    case 'Comsol'
        wd = pwd;
        cd(model.meshpath);
%        eval( ['cd ' model.meshpath ';' ]);
        eval( ['mesh = ' model.meshfile ';' ]);
        X = mesh.mesh('mesh1').getVertex';
        T = double( mesh.mesh('mesh1').getElem('tri')' +1 );
        mesh = TRI6( T, X, false );
        eval( ['cd ' wd ';' ]);
    
    % error
    otherwise
        error( 'unknown code' );
end
        