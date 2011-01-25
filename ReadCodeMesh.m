function mesh = ReadCodeMesh( model )
% READCODEMESH reads a mesh in a code-dependant format, and transforms it
% into a CArl format.
%
% syntax: mesh = ReadCodeMesh( model )

% R. Cottereau 01/2010

% switch on the type of code used
switch model.code
    
    % HOMEFE
    case 'HomeFE'
        mesh = model.mesh;

    % COMSOL
    case 'Comsol'
        mesh = struct( 'X', model.femcomsol.mesh.p', ...
                       'T', model.femcomsol.mesh.t(1:end-1,:)' );
    
    % error
    otherwise
        error( 'unknown code' );
end
        