function [Nm,Nc,opt] = Initiatevalues(model,coupling)
% function [Nm,Nc,mat,alpha1,alpha2,c2m,opt] = initialevalue(model,coupling)
% Set the initiale values for CArl and define opt (option) value.

% constants
Nm = length( model );    % number of models
Nc = length( coupling ); % number of coupling operations

% defaut values for option
opt.gerr = 1e-9 ; % error tolerance for geometry functions
opt.kappa = 1e-3 ; % kappa parameter for the coupling operator