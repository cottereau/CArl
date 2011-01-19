function M = MediatorSpace( mediator, Int, Rep )
% MEDIATORSPACE defines the weight matrix that projects the coupling
% matrices defined over the integration mesh onto the representation mesh
% chosen by the user
%
% syntax: M = MediatorSpace( mediator, Rep1, Rep2 )
%
%  mediator: structured array with the fields:
%       'type': 'deterministic' or 'stochastic' ('determinstic by default)
%       'support': 1 or 2. indicates which model should be used as the 
%                  support for the physical basis functions. Note that
%                  usually the coarser support should be used. In this 
%                  function only the field 'support' is used.
%  Rep: 1*2 cell containing the definition of the representation meshes of
%       the two models.
%
%  M: passage matrix from the integration mesh to the 
%              representation mesh 
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 05/2010

% choice of representation mesh
Rep = Rep{ mediator.support };

% constants
Ndof = size( Int.X, 1 );
Nrep = length( Rep.X );

% initialization
M = sparse( Ndof, Nrep );

% loop on elements of the representation mesh
for i1 = 1:Nrep
    M( Rep.Xr2Xi{i1}, i1 ) = Rep.value{i1};
end
