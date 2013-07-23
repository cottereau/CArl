function alpha = ArlequinWeight( cpl, i1, model )
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

% initializations
X = model.mesh.Points;

% define other element of the coupling
i2 = setdiff(1:2,i1);
val = cpl.cplval{i1};

% alpha in constant zones
ind1 = inside( cpl.free{i1}, X );
ind2 = inside( cpl.free12, X );
if model.mesh.d==1
    alpha = discontinuous1D( cpl.free{i1}, X(ind1,:), 1, ...
                             cpl.free12,   X(ind2,:), cpl.freeval{i1});
elseif model.mesh.d==2
    alpha = discontinuous( cpl.free{i1}, X(ind1,:), 1, ...
                           cpl.free12,   X(ind2,:), cpl.freeval{i1});
end

% alpha in the transition zone
% CONSTANT CASE
if numel(val)==1
    alpha = addRegion( alpha, cpl.domain, X, val );

% LINEAR CASE
elseif numel(val)==2
    ind = inside( cpl.domain, X );
    Xi = X(ind,:);
    d2 = cpl.free12.dist(Xi);
    d1 = cpl.free{i1}.dist(Xi);
    if all(isinf(d1)|isnan(d1))
        d1 = cpl.free{i2}.dist(Xi);
        f = (max(val)*d1 + min(val)*d2)./(d1+d2);
    elseif all(isinf(d2)|isnan(d2))
        d2 = cpl.free{i2}.dist(Xi);
        f = (max(val)*d2 + min(val)*d1)./(d1+d2);
    else
        f = (max(val)*d2 + min(val)*d1)./(d1+d2);
    end
    alpha = addRegion( alpha, cpl.domain, Xi, f );
end
   