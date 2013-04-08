classdef discontinuous
% Class for discontinous fonctions, defined as continuous functions over
% different subdomains
%
% vd = discontinuous(ls,v) defines a continuous function, defined as a 
% TriScatteredInterp class element over a single domain x defined as a 
% levelSet class element. Note that the nodes defining the function f do
% not have to be contained in the interior of the level-set. Equivalently,
% this means that the delaunay triangulation on which the interpolating
% function is defined is independent of the definition of the level-set.
%
% vd = discontinuous(ls,x,v) is equivalent to using 
% vd = discontinuous(ls,TriScatteredInterp(x,v)). Also you can set v as a
% single scalar.
% 
% vd = discontinuous(ls1,x1,v1,ls2,x2,v2,...,lsn,xn,vn) defines a 
% discontinuous function where the v1,v2,...,vn are defined over subdomains
% S1, S2, ..., SN. Note that all functions must be defined with 3 inputs.
%
% NB: no verification is made on the overlapping of subdomains
%
%  discontinuous properties:
%     x  - Ns*1 cell of subdomains (level-set class)
%     f  - Ns*1 cell of function definitions (TriScatteredInterp class)
%     Ns - number of subdomains
%
%  discontinuous methods:
%     interp    - Returns the values of the function at a series of points
%     plot      - plot the function over data set X
%     addRegion - plot the function over data set X
%
% see levelSet

% R. Cottereau 03/2013

    properties
        x    % Ns*1 cell of subdomains
        f    % Ns*1 cell of function definitions (TriScatteredInterp)
    end
    
    properties (Dependent)
        Ns;  % number of subdomains
    end
    
    methods
        % Definition of discontinuous function
        function obj = discontinuous(varargin)
            if nargin==2
                obj.x{1} = varargin{1};
                obj.f{1} = varargin{2};
            elseif nargin>=3
                obj.x{1} = varargin{1};
                if numel(varargin{3})==1
                    val = varargin{3}*ones(size(varargin{2},1),1);
                else
                    val = varargin{3};
                end
                obj.f{1} = TriScatteredInterp(varargin{2},val);
            end 
             if nargin>3 && mod(nargin,3)==0
                for i1 = 2:nargin/3
                    xi = varargin{(i1-1)*3+1};
                    fi = varargin{(i1-1)*3+2};
                    vali = varargin{(i1-1)*3+3};
                    obj = addRegion( obj, discontinuous( xi, fi, vali ) );
                end
             end
        end
        % setting of Ns property
        function N = get.Ns(obj)
            N = length(obj.x);
        end
        % plotting discontinuous function
        function plot( obj, X )
            [fval,ind] = interp( obj, X );
    %        fval = rand(size(fval))+fval;
            figure; hold on;
            for i1 = 1:obj.Ns
                if any(ind==i1)
                    dt = DelaunayTri(X(ind==i1,:));
                    xi = dt.X(:,1); yi = dt.X(:,2);
                    Xc = [mean(xi(dt.Triangulation),2) ...
                        mean(yi(dt.Triangulation),2)];
                    indT = distance( obj.x{i1}, Xc, true )<=0;
                    trisurf( dt.Triangulation(indT,:), dt.X(:,1), ...
                        dt.X(:,2), fval(ind==i1), fval(ind==i1) );
                    shading flat; view(2); colorbar;
                end
            end
        end
        % interpolation: getting values of f inside each subdomain
        function [fv,indg] = interp(obj,X)
            fv = zeros(size(X,1),1);
            indg = zeros(size(X,1),1);
            for i1 = 1:obj.Ns
                ind = distance(obj.x{i1},X,true)<=0;
                indg(ind) = i1;
                fv(ind) = obj.f{i1}(X(ind,:));
            end
        end
        % addregion: define a new subdomain and corresponding function
        function obj = addRegion(varargin)
            obj = varargin{1};
            if nargin==2
                obj1 = varargin{2};
                obj.x = [obj.x;obj1.x];
                obj.f = [obj.f;obj1.f];
            elseif nargin==4
                obj1 = discontinuous(varargin{2},varargin{3},varargin{4});
                obj = addRegion( obj, obj1 );
            end
        end
    end
end     