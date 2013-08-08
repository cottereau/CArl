function u = ArlequinOutput( model )
% ARLEQUINOUTPUT to construct solutions in format adapted to each of the
% models
%
%  u = ArlequinOutput( model )
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% brute solution
u = full(model.u);

% transformation depending on the model at hand
switch model.code
    
    % ACOUSTIC - HOMEFE
    case 'HomeFE'
        Ndof = size(model.HomeFE.mesh.X,1);
        u = mean( u( 1:Ndof, : ), 2 );
        
    % ACOUSTIC STOCHASTIC - MONTECARLOHOMEFE
    case 'MonteCarloHomeFE'
        Ndof = size(model.HomeFE.mesh.X,1);
        u = u( 1:Ndof, : );
        
    % ELASTIC - FE2D
    case 'FE2D'
        Ndof = size(model.FE2D.X,1);
        u = reshape( u( 1:2*Ndof, 1 ), 2, Ndof )';

    % COMSOL
    case 'Comsol'
        error('not implemented yet: output in COMSOL format')

    % TIMOSCHENKO BEAM - BEAM
    case 'Beam'
        Ndof = size(model.Beam.X,1);
        u = reshape( u( 1:2*Ndof, 1 ), Ndof, 2 );
                
    % unknown case
    otherwise
        error('this external code is not supported')
 
end

