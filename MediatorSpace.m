function M = MediatorSpace( mediator, Rep )
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

[ x, y, val ] = find( Rep{mediator}.M );
[~,~,x] = unique(x);
M = sparse( x, y, val, max(x), size(Rep{mediator}.M,2) );
