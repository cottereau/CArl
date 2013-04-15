function alpha = ArlequinWeight( cpl, i1, mesh )
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

% R. Cottereau 04/2010
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% define other element of the coupling
i2 = setdiff(1:2,i1);
val = cpl.cplval{i1};

% alpha in constant zones
alpha = discontinuous( cpl.free{i1}, mesh.X, 1, ...
                                      cpl.free12, mesh.X, cpl.freeval{i1});

% alpha in the transition zone
% CONSTANT CASE
if numel(cpl.cplval{i1})==1
    alpha = addRegion( alpha, cpl.mesh, mesh.X, cpl.cplval{i1} );

% LINEAR CASE
elseif numel(cpl.cplval{i1})==2
    ind = inside( cpl.mesh, mesh.X );
    Xi = mesh.X(ind,:);
    d2 = distance( cpl.free12, Xi, true );
    d1 = distance( cpl.free{i1}, Xi, true );
    if isempty(d1)
        d1 = distance( cpl.free{i2}, Xi, true );
        f = (max(val)*d1 + min(val)*d2)./(d1+d2);
    elseif isempty(d2)
        d2 = distance( cpl.free{i2}, Xi, true );
        f = (max(val)*d2 + min(val)*d1)./(d1+d2);
    else
        f = (max(val)*d2 + min(val)*d1)./(d1+d2);
    end
    alpha = addRegion( alpha, cpl.mesh, Xi, f );
end

