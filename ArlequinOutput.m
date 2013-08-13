function sol = ArlequinOutput( u, model, opt )
% ARLEQUINOUTPUT to construct solutions in format adapted to each of the
% models
%
%  u = ArlequinOutput( model )
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% initialization
Nm = length(model);
sol = cell(Nm,1);

% output the solution field in a format appropriate to each model
for i1 = 1:Nm
    
    % extract solution of the current model
    ui = full( u( (opt.iK(i1)+1):opt.iK(i1+1), : ) );

    % transformation depending on the model at hand
    switch model{i1}.code
        
        % ACOUSTIC - HOMEFE
        case 'HomeFE'
            Ndof = size(model{i1}.HomeFE.mesh.X,1);
            sol{i1} = mean( ui( 1:Ndof, : ), 2 );
            
            % ACOUSTIC STOCHASTIC - MONTECARLOHOMEFE
        case 'MonteCarloHomeFE'
            Ndof = size(model{i1}.HomeFE.mesh.X,1);
            sol{i1} = ui( 1:Ndof, : );
            
            % ELASTIC - FE2D
        case 'FE2D'
            Ndof = size(model{i1}.FE2D.X,1);
            sol{i1} = reshape( ui( 1:2*Ndof, 1 ), 2, Ndof )';
            
            % COMSOL
        case 'Comsol'
            error('not implemented yet: output in COMSOL format')
            
            % TIMOSCHENKO BEAM - BEAM
        case 'Beam'
            Ndof = size(model{i1}.Beam.X,1);
            sol{i1} = reshape( ui( 1:3*Ndof, 1 ), Ndof, 3 );
            
            % unknown case
        otherwise
            error('this external code is not supported')
            
    end
end

