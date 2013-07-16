function [ Nm, Nc, opt ] = Initiatevalues( model, coupling, opt )
% function [ Nm, Nc, opt ] = Initiatevalues( model, coupling, opt )
% Set the initiale values and default options.

% disabling some warnings
warning off MATLAB:DelaunayTri:ConsConsSplitWarnId;
warning off MATLAB:DelaunayTri:ConsSplitPtWarnId;
warning off MATLAB:TriRep:PtsNotInTriWarnId;
% constants
Nm = length( model );    % number of models
Nc = length( coupling ); % number of coupling operations

% error tolerance for geometry functions
if ~isfield( opt, 'gerr' )
    opt.gerr = 1e-9 ; 
end

% kappa parameter for the coupling operator
if ~isfield( opt, 'kappa' )
    opt.kappa = 1e-3;
end

% recompute coupling operator and matrices (the logicals correspond to each
% coupling)
if ~isfield( opt, 'recomputeC' )
    opt.recomputeC = true(Nc,1);
end

% recompute stiffness matrices (the logicals correspond to each model)
if ~isfield( opt, 'recomputeK' )
    opt.recomputeK = true(Nm,1);
end
