function [Nm,Nc,K,F,C1,C2,alpha1,alpha2,c2m,opt] = Initiatevalues(model,coupling)
% function [Nm,Nc,K,F,C1,C2,alpha1,alpha2,c2m,opt] = initialevalue(model,coupling)
% Set the initiale values for CArl and define opt (option) value.

% constants
Nm = length( model );    % number of models
Nc = length( coupling ); % number of coupling operations

% initialization
K = cell(Nm,1);
F = cell(Nm,1);
C1 = cell(Nc,1);
C2 = cell(Nc,1);
alpha1 = cell(Nc,1);
alpha2 = cell(Nc,1);
c2m = zeros(Nc,2);

% defaut values for option
opt.gerr = 1e-9 ; % error tolerance for geometry functions
opt.kappa = 1e-3 ; % kapap parameter for the coupling operator