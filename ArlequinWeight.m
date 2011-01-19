function [ alpha, distExt, distInt ] = ArlequinWeight( mesh, weight, n )
% ARLEQUINWEIGHT to compute the weight functions associated to the
% partition of models
%
% syntax: alpha = ArlequinWeight( model, coupling, n )
%
%  model: structured array containing the information relative to the 
%         model, in particular:
%       - 'type': describes the type of representation used for the
%                 geometry. This is used for all interpolation purposes,
%                 in particular for the definition of weight functions
%                 Implemented: {'FE' 'discrete'}
%       - 'mesh': array dependent on 'type'.
%  coupling: structured array describing the coupling options, in
%            particular:
%       - 'exteriorLevelSet': description of the outer boundary of the 
%                     coupling domain
%       - 'exteriorValue': 1*2 vector of the values of the weight function
%                     outside the exteriorLevelSet for each model. The sum
%                     of the values should be 1.
%       - 'interiorLevelSet': description of the inner boundary of the 
%                     coupling domain
%       - 'interiorValue': 1*2 vector of the values of the weight function
%                     inside the interiorLevelSet for each model. The sum
%                     of the values should be 1.
%       - 'weight': 'constant', 'linear'
%  n: indicates what model is being considered [1 or 2]
%
%  alpha: value of the weight function with the following format
%       - model.type='FE' vector of the same size as the number of lines in
%                    model.mesh.X
%       - model.type='discrete' not implemented yet
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2010

% constants
int = weight.interior;
ext = weight.exterior;

% test on the type of geometrical representation
switch mesh.type
    
    % FE representation
    case 'FE'
        
        % check coherency of weights
        int.value = [ int.value 1-int.value ];
        ext.value = [ ext.value 1-ext.value ];
        
        % intialization
        alpha = zeros( size(mesh.X,1), 1 );
        
        % computation of distances to interior and exterior
        distExt = DistanceMesh2LevelSet( mesh.X, ext.levelSet );
        distInt = - DistanceMesh2LevelSet( mesh.X, int.levelSet );
           
        % indices for the coupling domain
        ind = (distExt.*distInt) >= 0;

        % choice of the weight function
        switch weight.type
            
            
            % constant weight function
            case 'constant'
                alpha( ind ) = 1/2;
                
            % linear weight function
            case 'linear'
                a = abs(distInt(ind)) ./ (abs(distInt(ind))+abs(distExt(ind)));
                bInt = int.value(n);
                bExt = ext.value(n);
                alpha( ind ) = bInt + a*(bExt-bInt);
        
            % unknown weight function
            otherwise
                error( 'unknown type of weight function' )
        end

        % weight functions outside the coupling domain
        ind = ( distExt > 0 );
        alpha( ind ) = ext.value(n);
        ind = ( distInt > 0 );
        alpha( ind ) = int.value(n);

    % discrete representation
    case 'discrete'
        
        error( [ 'discrete representation for weight computation not ' ...
                 'implemented yet' ] );
             
    % unknown representation
    otherwise
        
        error( [ 'this type of representation for weight computation ' ...
                 'has not been implemented yet' ] );
end

