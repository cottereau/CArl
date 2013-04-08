classdef levelSet
% definition of a levelSet class for implicit definition of boundaries of
% complex convex or concave domains
%
% LS = levelSet(T,X) creates a level-set object based on the interface
% defined by the [Ne*2] connectivity matrix T and [Nn*2] coordinate matrix 
% X.
%
% LS = levelSet('circle',Xc,R) creates a level-set object based on a
% circular interface centered on Xc, with radius R.
%
% LS = levelSet('square',Xc,L) creates a level-set object based on a square
% interface with corner Xc and side length L.
%
% By default, the domain of interest is considered inside the interface
% controlled by the interface, but it is possible to consider the unbounded
% domain outside the interface by using LS = levelSet( ....,1)
%
%  levelSet properties:
%     X  - N*1 cell of Nn*2 coordinate matrix of nodes of the polygonal 
%          interfaces
%     T  - N*1 cell of Ne*2 connectivity matrices of polygonal interfaces
%     N  - number of disjoint level-sets
%     in - N*1 vector of logicals indicating the considered domains are 
%          inside (true) or outside (false) the different level-sets
%     sign - sign of the distance function inside the level-set
%
%  levelSet methods:
%     distance     - distance from interface to a set of points
%     interfaceI   - Returns a level-set made of only one closed surface
%     LSintersect  - intersecting two levelsets (new negative points are
%                    points that were formarly negative for both level-sets)
%     LScomplement - complement of two levelsets (new negative points are
%                    points that were in the first LS and not the second)
%     plot         - plot the distance function and interface
%     inside       - check whether points are inside the level-set

% NB: the structure of the class should be changed based on operational
% functions Union and Intersection ... for now, it is not robust at all for
% computing distances

% R. Cottereau 04/2013

    properties
        X     % Nn*2 coordinate matrix of nodes of the polygonal interface
        T     % Ne*2 connectivity matrices of polygonal interfaces
        in    % logical indicating the the domain is inside (true) or outside
              % (false) the levelset
    end

    properties (Constant)
        gerr = 1e-9; % error used to compute separations of points
        discC = 1e3; % discretization for approximation of circle convexhull 
    end
    
    properties (Dependent)
        sign  % sign of the distance function inside the level-set
        N     % number of disjoint level-sets
    end 
    
    methods
        % Definition of level-set
        function obj = levelSet(varargin)
            if ~ischar(varargin{1})
                obj.T{1} = varargin{1};
                obj.X{1} = varargin{2};
                obj.in = true;
                if nargin>2
                    obj.in = logical(varargin{3});
                end
            % special case of ellipse/circle
            elseif strcmp(varargin{1},'circle')
                theta = linspace(0,2*pi,obj.discC)';
                Xc = varargin{2};
                R = varargin{3};
                Xc = [Xc(:,1)+R*cos(theta) Xc(:,2)+R*sin(theta)];
                inc = true;
                if nargin>3
                    inc = varargin{4};
                end
                obj = levelSet( [(1:obj.discC-1)' (2:obj.discC)'], Xc, inc );
            % special case of rectangle/square
            elseif strcmp(varargin{1},'square')||strcmp(varargin{1},'rectangle')
                L = varargin{3};
                if length(L)==1
                    L = [L L];
                end
                Xr = repmat(varargin{2},[4 1]) + [0 0;L(1) 0;L;0 L(2)];
                inp = true;
                if nargin>3
                    inp = varargin{4};
                end
                obj = levelSet( [1 2;2 3;3 4;4 1], Xr, inp );
                return
            else
                error('unknown arguments')
            end
            if ~obj.in(1)
                obj.T{1} = obj.T{1}(end:-1:1,[2 1]);
            end
        end
        % get the sign of the distance function inside the level-set
        function s = get.sign(obj)
            s = -1*obj.in+1*~obj.in;
        end
        % get the number of disjoint level-sets
        function n = get.N(obj)
            n = length(obj.X);
        end
        % Distance from a set of points to the level set
        function d = distance( obj, X, val )
            Nx = size(X,1);
            d = zeros(Nx,obj.N);
            if nargin<3
                val = false;
            end
            for i1 = 1:obj.N
                poli = obj.X{i1}( [obj.T{i1}(:,1); obj.T{i1}(end,2) ], : );
                d(:,i1) = dppoli( X, poli );
                % change the sign for nodes inside the polygon
                ind = inpolygon( X(:,1), X(:,2), obj.X{i1}(:,1), obj.X{i1}(:,2));
                d(ind,i1) = obj.sign(i1)*d(ind,i1);
                d(~ind,i1) = -obj.sign(i1)*d(~ind,i1);
            end
            d = max(d,[],2);
            if ~val
                d = TriScatteredInterp( X, d );
            end
        end
        % get interface number i
        function lsi = interfaceI( obj, i1 )
            lsi = levelSet( obj.T{i1}, obj.X{i1}, obj.in(i1) );
        end
        % add an interface without checking for intersections
        function obj = addInterface( obj, T, X, in )
            n = obj.N+1;
            obj.X{n} = X;
            obj.T{n} = T;
            obj.in(n) = in;
        end
        % intersecting two domains
        function obj = LSintersect( obj, obj1 )
            for i1 = 1:obj.N
                for i2 = 1:obj1.N
                    x1 = obj.X{i1};
                    x2 = round(obj1.X{i2}/obj.gerr)*obj.gerr;
                    t1 = obj.T{i1};
                    t2 = obj1.T{i2};
                    nx = size(x1,1);
                    [in1,on1] = inpolygon( x1(:,1), x1(:,2), x2(:,1),x2(:,2) );
                    [in2,on2] = inpolygon( x2(:,1), x2(:,2), x1(:,1),x1(:,2) );
                    % normal case with two intersections
                    if ~all(in1&~on1)&&~all(in2&~on2)
                        t = [t1;t2+nx];
                        x = [x1;x2];
                        dt = DelaunayTri(x,t);
                        xc = dt.X(:,1); xc = mean(xc(dt.Triangulation),2);
                        yc = dt.X(:,2); yc = mean(yc(dt.Triangulation),2);
                        ind1 = inpolygon( xc, yc, x1(:,1), x1(:,2) );
                        if ~obj.in(i1)
                            ind1 = ~ind1;
                        end
                        ind2 = inpolygon( xc, yc, x2(:,1), x2(:,2) );
                        if ~obj1.in(i2)
                            ind2 = ~ind2;
                        end
                        xi = dt.X(:,1); yi = dt.X(:,2);
                        area = polyarea( xi(dt.Triangulation)', ...
                                               yi(dt.Triangulation)')'>eps;
                        ind = ind1&ind2&area;
                        dt = TriRep( dt.Triangulation(ind,:), dt.X );
                        [ff,xf] = freeBoundary(dt);
                        obj.T{i1} = ff;
                        obj.X{i1} = xf;
                        obj.in(i1) = true;
                    elseif i2==obj1.N
                        obj = addInterface( obj, t2, x2, obj1.in(i2) );
                    end
                end
            end
        end
        % creating a levelset as a complement of one level-set with respect
        % to another (points in the first LS and not the second)
        function obj = LScomplement( obj, obj1 )
            empty = false(obj.N,1);
            for i1 = 1:obj.N
                for i2 = 1:obj1.N
                    x1 = obj.X{i1};
                    x2 = round(obj1.X{i2}/obj.gerr)*obj.gerr;
                    t1 = obj.T{i1};
                    t2 = obj1.T{i2};
                    nx = size(x1,1);
                    [in1,on1] = inpolygon( x1(:,1), x1(:,2), x2(:,1),x2(:,2) );
                    [in2,on2] = inpolygon( x2(:,1), x2(:,2), x1(:,1),x1(:,2) );
                    % normal case
                    if ~all(in1&~on1)&&~all(in2&~on2)
                        t = [t1;t2+nx];
                        x = [x1;x2];
                        dt = DelaunayTri(x,t);
                        xc = dt.X(:,1); xc = mean(xc(dt.Triangulation),2);
                        yc = dt.X(:,2); yc = mean(yc(dt.Triangulation),2);
                        ind1 = inpolygon( xc, yc, x1(:,1), x1(:,2) );
                        if ~obj.in(i1)
                            ind1 = ~ind1;
                        end
                        ind2 = inpolygon( xc, yc, x2(:,1), x2(:,2) );
                        if ~obj1.in(i2)
                            ind2 = ~ind2;
                        end
                        xi = dt.X(:,1); yi = dt.X(:,2);
                        area = polyarea( xi(dt.Triangulation)', ...
                                               yi(dt.Triangulation)')'>eps;
                        ind = ind1&~ind2&area;
                        if all(~ind) && i2==obj1.N
                            empty(i1) = true;
                            break
                        elseif any(ind)
                            dt = TriRep( dt.Triangulation(ind,:), dt.X );
                            [ff,xf] = freeBoundary(dt);
                            obj.T{i1} = ff;
                            obj.X{i1} = xf;
                            obj.in(i1) = true;
                        end
                    elseif i2==obj1.N
                        obj = addInterface( obj, t2, x2, obj1.in(i2) );
                    end
                end
            end
            obj.T = obj.T(~empty);
            obj.X = obj.X(~empty);
            obj.in = obj.in(~empty);
        end
        % check whether points are inside the level-set
        function l = inside( obj, X )
            l = distance( obj, X, true )<=obj.gerr;
        end
        % plot function
        function  plot( obj, xv, yv )
            [x,y] = ndgrid(xv,yv);
            d = distance( obj, [x(:) y(:)], true );
            if ~isempty(d)
                figure; surf(xv,yv,reshape(d,length(xv),length(yv))')
                shading flat; view(2); colorbar
                hold on;
                contour( xv, yv, reshape(d,length(xv),length(yv))',[0 0], ...
                    'color','k','linewidth',2);
            else
                disp('empty level-set')
            end
        end
        function [n,c] = normal( varargin )
            obj = varargin{1};
            n = cell(obj.N,1);
            c = cell(obj.N,1);
            for i1 = 1:obj.N
                if nargin<2
                    ind = true(size(obj.T{i1},1),1);
                else
                    ind = varargin{2};
                end
                x = obj.X{i1}(:,1);
                y = obj.X{i1}(:,2);
                dx = x(obj.T{i1}(ind,2))-x(obj.T{i1}(ind,1));
                dy = y(obj.T{i1}(ind,2))-y(obj.T{i1}(ind,1));
                nn = [dy -dx];
                nn = sqrt(sum(nn.^2,2));
                n{i1} = [dy./nn -dx./nn];
                c{i1} = -n{i1}(:,1).*x(obj.T{i1}(ind,1)) ...
                        -n{i1}(:,2).*y(obj.T{i1}(ind,1));
            end
        end
    end
end
