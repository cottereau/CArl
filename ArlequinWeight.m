function alpha = ArlequinWeight( cpl, mesh1, mesh2 )
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

% alpha in constant zones
alpha{1} = discontinuous( cpl.free{1}, mesh1.X, 1, ...
                          cpl.free12, mesh1.X, cpl.freeval{1});
alpha{2} = discontinuous( cpl.free{2}, mesh2.X, 1, ...
                          cpl.free12, mesh2.X, cpl.freeval{2});

% alpha in the transition zone
if numel(cpl.cplval{1})==1      % constant case
    alpha{1} = addRegion( alpha{1}, cpl.mesh, mesh1.X, cpl.cplval{1} );
elseif numel(cpl.cplval{1})==2  % linear case
    ind = inside( cpl.mesh, mesh1.X );
    Xi = mesh1.X(ind,:);
    d2 = distance( cpl.free12, Xi, true );
    d1 = distance( cpl.free{1}, Xi, true );
    if isempty(d1)
        d1 = distance( cpl.free{1}, Xi, true );
        f = max(cpl.cplval{2})*d1./(d1+d2);
    else
        f = max(cpl.cplval{1})*d2./(d1+d2) + min(cpl.cplval{1})*d1./(d1+d2);
    end
    alpha{1} = addRegion( alpha{1}, cpl.mesh, Xi, f );
end
if numel(cpl.cplval{2})==1  % constant case
    alpha{2} = addRegion( alpha{1}, cpl.mesh, mesh2.X, cpl.cplval{2} );
elseif numel(cpl.cplval{2})==2  % linear case
    ind = inside( cpl.mesh, mesh2.X );
    Xi = mesh2.X(ind,:);
    d2 = distance( cpl.free12, Xi, true );
    d1 = distance( cpl.free{2}, Xi, true );
    if isempty(d1)
        d1 = distance( cpl.free{1}, Xi, true );
        f = max(cpl.cplval{2})*d1./(d1+d2);
    else
        f = max(cpl.cplval{1})*d2./(d1+d2) + min(cpl.cplval{1})*d1./(d1+d2);
    end
    alpha{2} = addRegion( alpha{2}, cpl.mesh, Xi, f );
end 
