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
    
    % raw output
    sol{i1} = model{i1}.u;
    
    % reconstruction of alpha.u
%     u = TriScatteredInterp( model{i1}.mesh.X3, model{i1}.u );
%     u = discontinuous( freeBoundary( model{i1}.mesh ), u );
%     model{i1}.alphaU = model{i1}.alpha.*u;
    
    % additional reconstruction for Monte Carlo solutions
    if isfield( model{i1}, 'uMC' );
        Nmc = size(model{i1}.uMC,2);
        model{i1}.auMC = zeros( [length(alpha) Nmc] );
        for i2 = 1:Nmc
            model{i1}.auMC(:,i2) = alpha .* model{i1}.uMC(:,i2);
        end
    end
end

