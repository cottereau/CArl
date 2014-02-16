function [ Nm, Nc, opt ] = Initiatevalues( model, coupling, opt )
% function [ Nm, Nc, opt ] = Initiatevalues( model, coupling, opt )
% Set the initial values and default options.

% disabling some warnings
warning off MATLAB:triangulation:PtsNotInTriWarnId
warning off MATLAB:delaunayTriangulation:ConsSplitPtWarnId
warning off MATLAB:delaunayTriangulation:ConsConsSplitWarnId
warning off MATLAB:delaunayTriangulation:DupPtsConsUpdatedWarnId
warning off MATLAB:delaunayTriangulation:DupConsWarnId
warning off MATLAB:scatteredInterpolant:DupPtsAvValuesWarnId

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

% compute solution
if ~isfield( opt, 'computeSol' )
    opt.computeSol = true;
end
