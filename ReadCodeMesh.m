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
    case {'HomeFE','MonteCarloHomeFE'}
        if ~isfield(model,'tri3')
            X = model.mesh.X;
            T = model.mesh.T;
            model.tri3 = struct('X',X,'Triangulation',T);
        else
            X = model.tri3.X;
            T = model.tri3.T;
        end
        d = size( X, 2 );
        if d==1
            model.mesh = struct( 'Triangulation', model.mesh.T, 'X', X );
        elseif d==2
            model.mesh = TRI6( T, X );
        else
            error('not implemented yet')
        end

    % COMSOL
    case 'Comsol'
        wd = pwd;
        cd(model.meshpath);
%        eval( ['cd ' model.meshpath ';' ]);
        eval( ['mesh = ' model.meshfile ';' ]);
        X = mesh.mesh('mesh1').getVertex';
        T = double( mesh.mesh('mesh1').getElem('tri')' +1 );
        mesh = TRI6( T, X );
        eval( ['cd ' wd ';' ]);
    
    % error
    otherwise
        error( 'unknown code' );
end
        