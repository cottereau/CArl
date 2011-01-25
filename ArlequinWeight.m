function alpha = ArlequinWeight( mesh, weight, LSet )
% ARLEQUINWEIGHT to compute the weight functions associated to the
% partition of models
%
% syntax: alpha = ArlequinWeight( mesh, weight, LSet )
%
%  mesh: structured array containing the fields 'X' and 'T'
%  weight: structured array containing the fields
%          -'value': value of the weight function between LSet.int and 
%                    LSet.ext, given as a vector of coefficients of a
%                    polynomial (a scalar for a constant weight, a 2*1
%                    vector for a linear weight, etc ...)
%          -'extvalue': value of the weight outside LSet.ext
%          -'intvalue': value of the weight inside the interior curve
%  n: indicates what model is being considered [1 or 2]
%
%  alpha: value of the weight function, given as a polynomial on each
%         element of the mesh [same size as T]
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2010

% warning
disp( [ 'for now, it is assumed that the boundary of the coupling zone' ...
        ' does not cross any element of any mesh' ] );

% intialization
alpha = zeros( size(mesh.T,1), size(weight.value,2) );

% value of level set function product at all nodes of elements

% alpha in the coupling domain
prodLS = LSet.ext.*LSet.int;
ind = any( prodLS( mesh.T )>= 1e-9, 2) | all( prodLS(mesh.T)==0, 2);
alpha(ind,:) = repmat( weight.value, [sum(ind) 1] );

% weight functions outside the coupling domain
ind = any( LSet.ext( mesh.T )>= 1e-9, 2);
alpha( ind, end ) = weight.extvalue;
ind = any( LSet.int( mesh.T )>= 1e-9, 2);
alpha( ind, end ) = weight.intvalue;

%         % linear weight function
%         a = abs(LSet.ext(ind)) ./ (abs(LSet.int(ind))+abs(LSet.ext(ind)));
%         bInt = weight.intvalue;
%         bExt = weight.extvalue;
%         alpha( ind ) = bInt + a*(bExt-bInt);
