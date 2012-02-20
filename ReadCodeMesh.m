function mesh = ReadCodeMesh( model )
% READCODEMESH reads a mesh in a code-dependant format, and transforms it
% into a CArl format.
%
% syntax: mesh = ReadCodeMesh( model )

% R. Cottereau 01/2010

% switch on the type of code used
switch model.code
    
    % HOMEFE
    case {'HomeFE','MonteCarloHomeFE'}
        X = model.mesh.X;
        d = size( X, 2 );
        if d==1
            mesh = struct( 'Triangulation', model.mesh.T, 'X', X );
        elseif d==2
            mesh = TriRep( model.mesh.T, X(:,1), X(:,2) );
        else
            error('not implemented yet')
        end

    % COMSOL
    case 'Comsol'
        mesh = struct( 'X', model.femcomsol.mesh.p', ...
                       'T', model.femcomsol.mesh.t(1:end-1,:)' );
    
    % error
    otherwise
        error( 'unknown code' );
end
        