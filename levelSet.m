classdef levelSet
% definition of a levelSet class for implicit definition of boundaries of
% complex convex or concave domains
%
% LS = levelSet(T,X) creates a level-set object considering a domain inside
% the interface defined by the [Ne*2] connectivity matrix T and [Nn*2] 
% coordinate matrix X.
%
% LS = levelSet('circle',Xc,R) creates a level-set object considering a
% domain inside a circular interface centered on Xc, with radius R.
%
% LS = levelSet('square',Xc,L) creates a level-set object considering a 
% domain inside a square interface with corner Xc and side length L.
%
% LS = levelSet(...,in) considers the interior of the level-set if in==true
% (default) and the outside of the level-set if in==false.
%
% LS = levelSet(true) corresponds to an empty domain and 
% LS = levelSet(false) corresponds to the full space.
%
%  levelSet properties:
%     X  - N*1 cell of Nn*2 coordinate matrix of nodes of the polygonal 
%          interfaces
%     T  - N*1 cell of Ne*2 connectivity matrices of polygonal interfaces
%     in - N*1 vector of logicals indicating the considered domains are 
%          inside (true) or outside (false) the different level-sets
%     N  - number of disjoint level-sets
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
%     normal       - normal vectors to edges in the level-set

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
            if nargin==1
                obj.T = cell(0);
                obj.X = cell(0);
                obj.in = varargin{1};
            elseif ~ischar(varargin{1})
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
        function lsi = getInterface( obj, i1 )
            lsi = levelSet( obj.T{i1}, obj.X{i1}, obj.in(i1) );
        end
        % replace the interface number i by another one
        function obj = setInterface( obj, obji, i1 )
            obj.X{i1} = obji.X{1};
            obj.T{i1} = obji.T{1};
            obj.in(i1) = obji.in(1);
        end
        % add an interface without checking for intersections
        function obj = addInterface( varargin )
            obj = varargin{1};
            if nargin==4
                n = obj.N+1;
                obj.X{n} = varargin{3};
                obj.T{n} = varargin{2};
                obj.in(n) = varargin{4};
            elseif nargin==2
                obj1 = varargin{2};
                obj = addInterface( obj, obj1.T{1}, obj1.X{1}, obj1.in );
            else
                error('leveSet/addInterface: incorrect number of arguments');
            end
        end
        % intersecting two domains composed of several level-sets
        function obj = intersection( obj1, obj2 )
            if (obj1.N==0&&obj1.in)||(obj2.N==0&&obj2.in)
                obj = levelSet(true);
                return
            end
            if (obj1.N==0&&~obj1.in)
                obj = obj2;
                return
            end
            if (obj2.N==0&&~obj2.in)
                obj = obj1;
                return
            end
            ls1 = getInterface( obj1, 1 );
            ls2 = getInterface( obj2, 1 );
            obj = LSintersection( ls1, ls2 );
            if obj1.N>1||obj2.N>1
                for i1 = 1:obj1.N
                    for i2 = 1:obj2.N
                        ls1 = getInterface( obj1, i1 );
                        ls2 = getInterface( obj2, i2 );
                        obj = addInterface( obj, LSintersection( ls2, ls1 ));
                    end
                end
                obj = collapseIntersection( obj );
            end
        end
        % union of two domains composed of several level-sets
        function obj = union( obj1, obj2 )
            if (obj1.N==0&&~obj1.in)||(obj2.N==0&&~obj2.in)
                obj = levelSet(false);
                return
            end
            if (obj1.N==0&&obj1.in)
                obj = obj2;
                return
            end
            if (obj2.N==0&&obj2.in)
                obj = obj1;
                return
            end
            ls1 = getInterface( obj1, 1 );
            ls2 = getInterface( obj2, 1 );
            obj = LSunion( ls1, ls2 );
            if obj1.N>1||obj2.N>1
                for i1 = 1:obj1.N
                    for i2 = 1:obj2.N
                        ls1 = getInterface( obj1, i1 );
                        ls2 = getInterface( obj2, i2 );
                        obj = addInterface( obj, LSunion( ls2, ls1 ));
                    end
                end
            end
        end
        % returns the intersection of all the level-sets in one object
        function obj = collapseIntersection(obj)
            for i1 = obj.N:-1:2
                for i2 = obj.N-1:-1:1
                    ls1 = getInterface( obj, i1 );
                    ls2 = getInterface( obj, i2 );
                    ls12 = LSintersection( ls2, ls1 );
                    if ls12.N<2
                        obj = setInterface( obj, ls12, i1 );
                    end
                end
            end
        end
        % returns the union of all the level-sets in one object
        function obj = collapseUnion(obj)
            for i1 = obj.N:-1:2
                for i2 = obj.N-1:-1:1
                    ls1 = getInterface( obj, i1 );
                    ls2 = getInterface( obj, i2 );
                    ls12 = LSunion( ls2, ls1 );
                    if ls12.N<2
                        obj = setInterface( obj, ls12, i1 );
                    end
                end
            end
        end
        % intersecting two level-sets of size 1 each
        function obj = LSintersection( obj1, obj2 )
            l12 = inside( obj1, obj2.X{1}, false );
            l21 = inside( obj2, obj1.X{1}, false );
            if ~all(l12)&&~all(l21)
                [dt,ind1,ind2] = CrossLS( obj1, obj2 );
                dt = subSet( dt, ind1&ind2 );
                obj = freeBoundary(dt);
            elseif all(l12)&&all(l21)
                obj = addInterface( obj1, obj2 );
            elseif all(l21)
                obj = obj1;
            elseif all(l12)
                obj = obj2;
            else
                obj = levelSet(true);
            end
        end
        % union of two level-sets of size 1 each
        function obj = LSunion( obj1, obj2 )
            l12 = inside( obj1, obj2.X{1}, false );
            l21 = inside( obj2, obj1.X{1}, false );
            if ~all(l12)&&~all(l21)
                [dt,ind1,ind2] = CrossLS( obj1, obj2 );
                dt = subSet( dt, ind1|ind2 );
                obj = freeBoundary(dt);
            elseif all(l21)&&~all(l12)
                obj = obj2;
            elseif all(l12)&&~all(l21)
                obj = obj1;
            else
                obj = levelSet(false);
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
        function l = inside( obj, X, lon )
            if nargin<3 || lon
                l = distance( obj, X, true )<=obj.gerr;
            else
                l = distance( obj, X, true )<=-obj.gerr;
            end
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
        % normal vectors to edges in the level-set
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
        % create a Delaunay triangulation constrained by two level-sets
        % (always convex) and indicates which elements are in/out 
        function [dt,ind1,ind2] = CrossLS( obj1, obj2, i1, i2 )
            if nargin<4 || isempty(i2);
                i2 = 1;
            end
            if nargin<3 || isempty(i1);
                i1 = 1;
            end
            C = [ obj1.T{i1}; obj2.T{i2}+size(obj1.X{i1},1) ];
            Xdt = [ obj1.X{i1}; obj2.X{i2} ];
            dt = DelaunayTri( Xdt, C );
            dt = TRI6( dt.Triangulation, dt.X );
            ind1 = inside( obj1, dt.X, true );
            ind1 = all( ind1(dt.T), 2 );
            ind2 = inside( obj2, dt.X, true );
            ind2 = all( ind2(dt.T), 2 );
        end
        % get rid of nodes that are not used in T, and of repeated elements
        function obj = cleanT(obj)
            % get rid of elements that are repeated
            [~,ind] = unique( sort(obj.T{1},2) ,'rows' );
            obj.T{1} = obj.T{1}(ind,:);
            % get rid of elements that join the same nodes
            ind = diff(obj.T{1},[],2)==0;
            obj.T{1} = obj.T{1}(~ind,:);
            % get rid of nodes that are not used in T
            ind = unique(obj.T{1});
            n3 = length(ind);
            for i1 = 1:n3
                obj.T{1}( obj.T{1}==ind(i1) ) = i1;
            end
            obj.X{1} = obj.X{1}(ind,:);
        end
        % get rid of nodes that are repeated
        function obj = cleanX(obj)
            xrnd = round(obj.X{1}/obj.gerr)*obj.gerr;
            [~,indx,indu] = unique( xrnd, 'rows' );
%            [~,indx,indu] = unique( xrnd, 'rows', 'stable' );
            Nx = size(obj.X{1},1);
            if length(indx)<Nx
                obj.X{1} = obj.X{1}(indx,:);
                obj.T{1} = indu(obj.T{1});
            end
        end
        % clean X of unused nodes and repeated nodes and elements
        function obj = clean(obj)
            for i1 = obj.N
                obji = interfaceI( obj, i1 );
                obji = cleanX(obji);
                obji = cleanT(obji);
                obj.X{i1} = obji.X{1};
                obj.T{i1} = obji.T{1};
            end
        end
    end
end
