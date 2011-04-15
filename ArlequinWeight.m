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
%                    For linear : give the value as :
%              [related to interior domain , related to exterior domain]
%          -'extvalue': value of the weight outside LSet.ext
%          -'intvalue': value of the weight inside the interior curve
%
%  alpha: value of the weight function, given as a polynomial on each
%         element of the mesh [same size as T]
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% R. Cottereau 04/2010
% C. Zaccardi 04/2011
        
% intialization
alpha = zeros( size(mesh.T) );

% weight functions outside the coupling domain (outside exterior LSet)
ind = any( LSet.ext( mesh.T )<= -1e-9, 2);
alpha( ind,: ) = weight.intvalue ;
ind = any( LSet.ext( mesh.T )>= 1e-9, 2);
alpha( ind,: ) = weight.extvalue ;

% weight functions outside the coupling domain (inside exterior LSet)
ind = any( LSet.int( mesh.T )<= -1e-9, 2);
alpha( ind,: ) = weight.extvalue ;
ind = any( LSet.int( mesh.T )>= 1e-9, 2);
alpha( ind,: ) = weight.intvalue ;

% alpha in the coupling domain
prodLS = LSet.ext.*LSet.int;
indT = any( prodLS( mesh.T )>= 1e-9, 2) | all( prodLS(mesh.T)==0, 2);
[indX,j1,j2] = unique(mesh.T(indT,:));
a = abs(LSet.int(indX)) ./ (abs(LSet.ext(indX))+abs(LSet.int(indX)));
if size(weight.value,2)==2
    a = weight.value(1) + a*(weight.value(2)-weight.value(1));
else
    error('Not implemented yet')
end
alpha( indT,: ) = reshape( a(j2), sum(indT), size(mesh.T,2) );





