function [ sol, model ] = ArlequinOutput( model )
% ARLEQUINOUTPUT to construct an Arlequin solution based on the individual
% solutions for each model, and the alpha functions
%
%  u = alpha*u

% developed at 
% Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579, 
% grande voie des vignes
% F-92295 Chatenay-Malabry
% FRANCE
% contact: regis.cottereau@ecp.fr

% constants
Nm = length(model);

% initialization
sol = cell(Nm,1);

% construction of local solutions to each model
for i1 = 1:Nm
    
    % solution brute
    sol{i1} = full(model{i1}.u);
    
    % additional reconstruction for Monte Carlo solutions
    if isfield( model{i1}, 'uMC' );
        sol{i1} = mean( full(model{i1}.uMC), 2 );
    end
    
    % reconstruction of alpha.u
    alpha = model{i1}.alpha;
    u = TriScatteredInterp( model{i1}.mesh.X3, sol{i1} );
    u = discontinuous( freeBoundary(model{i1}.mesh), u );
%    model{i1}.au = alpha.*u;
    
end

