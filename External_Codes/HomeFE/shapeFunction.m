function [N,Nxi,Neta] = shapeFunction( elem, nen, pospg ) 
% [N,Nxi,Neta] = shapeFunction( elem, nen, pospg ) 
%
% Evaluate the shape functions in the gauss points 
%
% INPUT
%   elem          type of element (0: interval, 1: triangle, 2:square)
%   nen           number of nodes per element
%   pospg         position of the gauss points in the reference
%
% OUTPUT
%   N,Nxi,Neta    matrices with the shape functions and their derivatives
%                 evaluated in the gauss points. For the interval, Neta=[]
%
if elem == 0 % interval
    xi = pospg; 
    if nen == 1  %P0
      N = ones( size( xi ) );
      Nxi = zeros( size( xi ) );
    elseif nen == 2  %P1
      N = [ 1-xi, xi ]; 
      Nxi = [ -ones(size(xi)) ones(size(xi)) ]; 
    else
      error( 'shapeFunction error: unknown number of nodes for interval' )
    end
    Neta = [];
        
elseif elem == 1  % triangles
    xi = pospg(:,1); 
    eta = pospg(:,2); 
    if nen == 1  %P0
      N = ones( size( xi ) );
      Nxi = zeros( size( xi ) );
      Neta = zeros( size( xi ) );
    elseif nen == 3  %P1
      N = [xi,eta,1-(xi+eta)]; 
      Nxi = [ones( size( xi ) ), zeros( size( xi ) ),-ones( size( xi ) )]; 
      Neta = [zeros( size( xi ) ), ones( size( xi ) ),-ones( size( xi ) )]; 
    elseif nen == 6  %P2
      N = [xi.*(2*xi-1),eta.*(2*eta-1),(1-2*(xi+eta)).*(1-(xi+eta)), ...
           4*xi.*eta,4*eta.*(1-(xi+eta)),4*xi.*(1-(xi+eta))]; 
      Nxi = [4*xi-1,zeros(size(xi)),-3+4*(xi+eta),4*eta, ...
            -4*eta,4*(1-2*xi-eta)]; 
      Neta = [zeros(size(xi)),4*eta-1,-3+4*(xi+eta),4*xi, ...
              4*(1-xi-2*eta),-4*xi]; 
    elseif nen == 4  %P1+
      N = [xi,eta,1-(xi+eta),27*xi.*eta.*(1-xi-eta)]; 
      Nxi = [ones( size( xi ) ), zeros( size( xi ) ), ...
         -ones( size( xi ) ), 27*eta.*(1-2*xi-eta)]; 
      Neta = [zeros( size( xi ) ), ones( size( xi ) ), ...
         -ones( size( xi ) ), 27*xi.*(1-2*eta-xi)]; 
    else
      error( 'shapeFunction error: triangle not implemented' )
    end
   
elseif elem == 2   % squares
    xi = pospg(:,1); 
    eta = pospg(:,2); 
    if nen == 1  % P0
      N = ones( size( xi ) );
      Nxi = zeros( size( xi ) );
      Neta = zeros( size( xi ) );
    elseif nen == 4  % Q1
      N = [(1-xi).*(1-eta)/4, (1+xi).*(1-eta)/4, ...
           (1+xi).*(1+eta)/4, (1-xi).*(1+eta)/4]; 
      Nxi = [(eta-1)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4]; 
      Neta = [(xi-1)/4, -(1+xi)/4,   (1+xi)/4,  (1-xi)/4 ]; 
    elseif nen == 9  %Q2
         N = [xi.*(xi-1).*eta.*(eta-1)/4, xi.*(xi+1).*eta.*(eta-1)/4, ...
              xi.*(xi+1).*eta.*(eta+1)/4, xi.*(xi-1).*eta.*(eta+1)/4, ...
              (1-xi.^2).*eta.*(eta-1)/2,  xi.*(xi+1).*(1-eta.^2)/2,   ...
              (1-xi.^2).*eta.*(eta+1)/2,  xi.*(xi-1).*(1-eta.^2)/2,   ...
              (1-xi.^2).*(1-eta.^2)];
         Nxi = [(xi-1/2).*eta.*(eta-1)/2,   (xi+1/2).*eta.*(eta-1)/2, ...
                (xi+1/2).*eta.*(eta+1)/2,   (xi-1/2).*eta.*(eta+1)/2, ...
                -xi.*eta.*(eta-1),          (xi+1/2).*(1-eta.^2),   ...
                -xi.*eta.*(eta+1),          (xi-1/2).*(1-eta.^2),   ...
                -2*xi.*(1-eta.^2)];
         Neta = [xi.*(xi-1).*(eta-1/2)/2,    xi.*(xi+1).*(eta-1/2)/2, ...
                 xi.*(xi+1).*(eta+1/2)/2,    xi.*(xi-1).*(eta+1/2)/2, ...
                 (1-xi.^2).*(eta-1/2),       xi.*(xi+1).*(-eta),   ...
                 (1-xi.^2).*(eta+1/2),       xi.*(xi-1).*(-eta),   ...
                 (1-xi.^2).*(-2*eta)];
    else
        error( 'shapeFunction error: square not implemented' )
    end
else 
   error( 'shapeFunction error' )
end

