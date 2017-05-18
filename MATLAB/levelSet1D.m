classdef levelSet1D
% definition of a levelSet class for implicit definition of boundaries of
% complex convex or concave domains in 1D
%
% LS = levelSet1D(x1,x2) creates a 1D domain between x1 and x2
%
% LS = levelSet1D(...,in) considers the interior of [x1 x1] if in==true
% (default) and the outside if in==false.
%
% LS = levelSet1D(true) corresponds to an empty domain and 
% LS = levelSet1D(false) corresponds to the full line.
%
% When there are several segments defined, the composite domain is defined
% as the union of all the segments
%
%  levelSet properties:
%     X  - N*2 matrix of coordinates of nodes of each of the domains
%     N  - number of disjoint domains
%     sign - sign of the distance function inside the level-set
%
%  levelSet methods:
%     getInterface - Returns a single interval subdomain 
%     setInterface - set a particular interval in a composite domain
%     rmInterface  - remove one interval from a composite domain
%     intersection - intersection of two domains
%     union        - union of two domains
%     complement   - complement of one domain within another
%     distance     - distance from a set of points to the closest boundary
%     boundary     - find the point on the boundary of the domain
%     inside       - check whether points are inside the domain

% R. Cottereau 07/2013

    properties
        X     % N*2 coordinate matrix of nodes of the domains
        dir   % direction of the 1D levelSet embedded in a 2D space
        x0    % reference point of the 1D levelSet embedded in 2D space
    end

    properties (Constant)
        gerr = 1e-9; % error used to compute separations of points
        sign = -1;   % value of the distance function inside the domain
    end
    
    properties (Dependent)
        N     % number of disjoint domains
    end 
    
    methods
        % Definition of level-set
        function obj = levelSet1D( x1, x2, in, dir, x0 )
            if nargin==1
                obj = levelSet1D( -Inf, Inf, x1 );
            elseif nargin==2
                obj = levelSet1D( x1, x2, true );
            else
                ind = x2>=x1;
                x = sort( [x1(ind,:) x2(ind,:)], 2, 'ascend' );
                obj.X = sortrows( x, 1 );
                if ~in
                    obj.X = [-Inf         x1(1)      ; 
                              x2(1:end-1) x1(2:end) ;
                              x2(end)     Inf        ];
                    if all(isinf(obj.X(1,:)))
                        obj = rmInterface( obj, 1 );
                    end
                    if all(isinf(obj.X(end,:)))
                        obj = rmInterface( obj, obj.N );
                    end
                end
            end
            if nargin<4 || isempty(dir)
                obj.dir = [1 0];
            else
                obj.dir = dir;
            end
            if nargin<5 || isempty(x0)
                obj.x0 = [0 0];
            else
                obj.x0 = x0;
            end            
            obj = clean(obj);
        end
        % get the number of disjoint level-sets
        function n = get.N(obj)
            n = size(obj.X,1);
        end
        % get interface number i
        function lsi = getInterface( obj, i1 )
            lsi = levelSet1D( obj.X(i1,1), obj.X(i1,2) );
        end
        % erase an interface of a levelSet structure
        function obj = rmInterface( obj, i1 )
            obj.X = obj.X([1:(i1-1) (i1+1):obj.N],:);
        end
        % replace the interface number i by another one
        function obj = setInterface( obj, obji, ind )
            if nargin<3
                ind = obj.N+(1:obji.N);
            end
            x = obj.X;
            x(ind,:) = obji.X;
            obj = levelSet1D( x(:,1), x(:,2) );
        end
        % complement of obj2 in obj1
        function obj = complement( obj1, obj2 )
            obj2 = levelSet1D( obj2.X(:,1), obj2.X(:,2), false );
            obj = intersection( obj1, obj2 );
        end
        % intersecting two domains
        function obj = intersection( obj1, obj2 )
            if obj1.N==1 && obj2.N==1
                x1 = max( obj1.X(1), obj2.X(1) );
                x2 = min( obj1.X(2), obj2.X(2) );
                if ~isempty(x1) && abs(x1-x2)>obj1.gerr
                    obj = levelSet1D( x1, x2, true );
                else
                    obj = levelSet1D( false );
                end
            else
                obj = levelSet1D( false );
                for i1 = 1:obj1.N
                    for i2 = 1:obj2.N
                        ls1 = getInterface( obj1, i1 );
                        ls2 = getInterface( obj2, i2 );
                        obj = union( obj, intersection( ls1, ls2 ) );
                    end
                end
            end
        end
        % union of two domains composed of several level-sets
        function obj = union( obj, obj2 )
            obj = setInterface( obj, obj2 );
            obj = clean( obj );
        end
        % Distance from a set of points to the level set
        function [d,ind] = dist( obj, X )
            if size(X,1)<size(X,2)
                X = X';
            end
            d = Inf(size(X));
            ind = zeros(size(X));
            for i1 = 1:obj.N
                d1 = abs( X - obj.X(i1,1) );
                d2 = abs( X - obj.X(i1,2) );
                d = min( [d d1 d2], [], 2 );
                ind( X>=obj.X(i1,1) & X<=obj.X(i1,2)) = i1;
            end
            lind = ind>0;
            d(~lind) = -obj.sign .* d(~lind);
            d(lind) = obj.sign .* d(lind);
            if any(isinf(d))
                d = [];
            end
        end
        % check whether points are inside the level-set
        function [in,on] = inside( obj, X, lon )
            d = dist( obj, X );
            on = abs(d)<= obj.gerr;
            in = d <= -obj.gerr;
            if nargin<3 || lon
                in = in|on;
            end
        end
        % find the point on the boundary of the domain
        function x = boundary( obj, x1, x2 )
            d1 = dist( obj, x1 );
            d2 = dist( obj, x2 );
            dx = abs(d1.*(x2-x1)./(d2-d1));
            ind = x1<x2;
            x(ind) = x1(ind) + dx(ind);
            x(~ind) = x1(~ind) + dx(~ind);
        end
        % remove unncessary segments (merge segments partially overlapping
        % and fully included in another one). remove empty subdomains
        function obj = clean( obj )
            N = obj.N;
            for i1 = N:-1:2
                if any(obj.X(1:i1-1,2)>=obj.X(i1,2))
                    obj = rmInterface( obj, i1 );
                elseif any(obj.X(1:i1-1,2)>=obj.X(i1,1))
                    obj.X(i1-1,2) = obj.X(i1,2);
                    obj = rmInterface( obj, i1 );
                end
            end
        end
    end
end
