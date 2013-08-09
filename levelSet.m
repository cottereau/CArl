classdef levelSet
% definition of a levelSet class for implicit definition of domains
%
% LS = levelSet(X,V) creates a level-set object, from the data of a set of
% nodes X and the value of the distance function V evaluated at these
% nodes.
%
% LS = levelSet(X,'circle',Xc,R) creates a level-set object considering a
% domain inside a circular interface centered on Xc, with radius R. The
% level set is evaluated at points x.
%
% LS = levelSet(X,'square',Xc,L) creates a level-set object considering a
% domain inside a square interface with corner Xc and side length L.
% LS = levelSet(X,'rectangle',Xc,L) is the same with L a vector of two
% coordinates
%
% LS = levelSet(X,'polygon',Xp,Tp) creates a level-set object 
% considering a polygon defined by nodes Xp, and connectivity Tp.
%
% LS = levelSet(X,'init',d,Tp) creates a level-set object 
% given a boundary defined by nodes X, and connectivity Tp, and distance 
% function d (the distance is not recomputed and no new points are added).
%
% LS = levelSet(...,in) considers the exterior of the interface if 
% in==false, and the interior if in==true (default).
%
% LS = levelSet(true) corresponds to an empty domain and
% LS = levelSet(false) corresponds to the full space.
% LS = levelSet is the same as levelSet(true)
%
%    levelSet properties:
%        mesh    - distance function from interface [scatteredInterpolant]
%        dist    - corresponding triangulation [delaunayTriangulation]
%        d       - dimension of space [scalar]
%
%    levelSet methods:
%     intersection - intersection of two domains
%     union        - union of two domains
%     complement   - complement of one domain within another
%     inside       - check whether points are inside the domain
    
    properties
        dist  % distance function from interface [scatteredInterpolant]
        mesh  % corresponding triangulation [delaunayTriangulation]
    end
    
    properties (Constant)
        gerr = 1e-9; % error used to compute separations of points
        in = -1;     % value of the distance function inside the domain
    end
    
    methods
        % Definition of level-set
        function obj = levelSet(varargin)
            X = [];
            V = [];
            C = [];
            if nargin==0
            elseif nargin==1
                X = [0 0;1 0;1 1];
                V = obj.in*inf(3,1);
            elseif nargin==2||nargin==3
                X = varargin{1};
                V = varargin{2};
            else
                x = varargin{1};
                geom = varargin{2};
                x1 = varargin{3};
                x2 = varargin{4};
                switch geom
                    % CIRCLE
                    case 'circle'
                        nc = 1e3;
                        t = linspace( 0, 2*pi, nc )';
                        xon = [x1(1)+x2*cos(t) x1(2)+x2*sin(t)];
                        xout = [x1(1)+2*x2*cos(t) x1(2)+2*x2*sin(t)];
                        if ~isempty(x)
                            v = x2 - sqrt((x(:,1)-x1(1)).^2+(x(:,2)-x1(2)).^2);
                        else
                            v = [];
                        end
                        X = [ xon; x; x1; xout ];
                        V = obj.in * [ zeros(nc,1); v; x2; -x2*ones(nc,1) ];
                        C = [(1:nc-1)' (2:nc)'];
                        
                    % RECTANGLE OR SQUARE
                    case {'square','rectangle'}
                        xm = min(x2)/2;
                        de = 0.9*xm;
                        if numel(x2)==1
                            x2 = [x2 x2];
                        end
                        xv = x1(1) + [-de 0 x2(1)/2 x2(1) x2(1)+de];
                        yv = x1(2) + [-de 0 x2(2)/2 x2(2) x2(2)+de];
                        [xv,yv] = ndgrid(xv,yv);
                        xp = [xv(:) yv(:)];
                        xc = x1(1) + [xm x2(1)-xm]';
                        yc = x1(2) + [xm x2(2)-xm]';
                        X = [ xp; [xc yc]; x ];
                        ind = [7 8 9 14 19 18 17 12];
                        C = [ind' [ind(2:end) ind(1)]'];
                        V = -obj.in * dppoli( X, C, X );
                        ind = inpolygon( X(:,1), X(:,2), xp(ind,1), xp(ind,2) );
                        V(ind) = -V(ind);
                        
                    % POLYGON
                    case 'polygon'
                        de = 0.4;
                        e1 = x1(x2(:,1),:);
                        e2 = x1(x2(:,2),:);
                        dir = e2 - e1;
                        dir = de * [dir(:,2) -dir(:,1)];
                        xv = x1(x2(:,1),1);
                        yv = x1(x2(:,1),2);
                        xsup = [e1+dir; e1-dir; e2+dir; e2-dir];
                        X = [ x1; x; xsup ];
                        V = -obj.in * dppoli( X, x2, x1 );
                        ind = inpolygon( X(:,1), X(:,2), xv, yv );
                        V(ind) = -V(ind);
                        C = x2;

                    % RE-INITIALIZATION
                    case 'init'
                        X = x;
                        C = x2;
                        V = x1;

                    otherwise
                        error('unknown type of geometrical domain')
                end
            end
            
            % create structure
            if ~isempty(V)
                obj.dist = scatteredInterpolant( X, V, 'linear', 'nearest' );
            else
                obj.dist = scatteredInterpolant;
            end
            if ~isempty(C)
                obj.mesh = delaunayTriangulation( X, C );
                obj.mesh = delaunayTriangulation( obj.mesh.Points, ...
                                                    obj.mesh.Constraints );
            else
                obj.mesh = delaunayTriangulation(X);
            end
            
            % complement
            if (nargin==1||nargin==3||nargin==5) && ~varargin{end}
               obj = complement(obj);
            end
        end
        % plot function
        function  plot(obj)
            if all(isinf(obj.dist.Values))
                if all(obj.in*obj.dist.Values>0)
                    disp('full domain')
                elseif all(obj.in*obj.dist.Values<0)
                    disp('empty domain')
                end
                return
            end
            T = obj.mesh.ConnectivityList;
            X = obj.dist.Points;
            xc = incenter(obj.mesh);
            vc = obj.dist(xc);
            figure; trimesh( T, X(:,1), X(:,2), 'line','-','color','k' );
            hold on; scatter( xc(:,1), xc(:,2), 50, vc, 'full' )
            hold on; scatter( X(:,1), X(:,2), 50, obj.dist.Values )
            [ff,xf] = boundary(obj);
            x = [xf(ff(:,1),1) xf(ff(:,2),1)];
            y = [xf(ff(:,1),2) xf(ff(:,2),2)];
            hold on; plot( x', y', 'k-', 'LineWidth', 2 );
            colorbar
        end
        % get the interface (for which distance=0)
        function varargout = boundary(obj)
            C = obj.mesh.Constraints;
            % both nodes of the constraint are on the boundary
            [~,on] = inside( obj, obj.dist.Points(C(:),:), true );
            C = C( all( reshape(on,size(C,1),2), 2 ), : );
            % a node inside the domain is not on two segments
            sc = sort(C(:));
            dn = diff([ -Inf; sc; Inf ]);
            dn = dn(2:end)==0 | dn(1:end-1)==0;
            C = C( all( ismember(C,sc(dn)), 2 ), : );
            % attached nodes are on different sides of the interface (or
            % on the boundary of the convex hull)
            on = true(size(C,1),1);
            bnd = ~all( ismember(C,convexHull(obj.mesh)), 2 );
            if ~isempty(bnd)
                TE = edgeAttachments( obj.mesh, C(bnd,1), C(bnd,2) );
                TE = cat( 1, TE{:} )';
                XE = obj.mesh.ConnectivityList( TE(:), :);
                XE = reshape( XE', 6, size(XE,1)/2 )';
                VE = obj.dist.Values(XE);
                VE( abs(VE)<obj.gerr ) = 0;
                VE = sum( sign( VE ), 2 );
                on(bnd) = ~( VE==2 | VE==-2 );
            end
            C = C(on,:);
            % output
            varargout{1} = C;
            if nargout==2
                [xf,~,ff] = unique(varargout{1});
                varargout{1} = reshape(ff,size(varargout{1}));
                varargout{2} = obj.mesh.Points(xf,:);
            end
        end
        % complement of obj1 in obj
        function obj = complement( obj, obj1 )
            if nargin==1
                obj.dist.Values = -obj.dist.Values;
            else
                obj = intersection( obj, complement(obj1) );
            end
        end
        % intersecting two domains
        function obj = intersection(obj1,obj2)
            C1 = obj1.mesh.Constraints( reduceBoundary(obj1,obj2), : );
            C2 = obj2.mesh.Constraints( reduceBoundary(obj2,obj1), : );
            obj1.mesh = delaunayTriangulation( obj1.mesh.Points, C1 );
            obj2.mesh = delaunayTriangulation( obj2.mesh.Points, C2 );
            dt = mergeMeshes( obj1, obj2 );
            if isempty(C1)&&isempty(C2)
                V = -obj1.in*inf(size(dt.Points,1),1);
            else
                V = max([obj1.dist(dt.Points) obj2.dist(dt.Points)],[],2);
            end
            obj = levelSet( dt.Points, 'init', V, dt.Constraints );
            obj = clean(obj);
        end
        % check whether points are inside the domain
        function [in,on] = inside( obj, X, lon )
            d = obj.dist( X );
            on = abs(d)<= obj.gerr;
            in = d <= -obj.gerr;
            if nargin<3 || lon
                in = in|on;
            end
        end
        % clean level-Set
        function obj = clean(obj)
            if all( inside(obj,obj.mesh.Points,true) )
                obj = levelSet(true);
            elseif all( ~inside(obj,obj.mesh.Points,false) )
                obj = levelSet(false);
            else
                obj = updateDistance(obj);
                obj = updateDistance(obj);
            end
        end
        % update distance function in level-set
        function obj = updateDistance(obj)
            C = boundary(obj);
            X = obj.mesh.Points;
            d = dppoli( X, C, X );
            [inin,on] = inside( obj, X, false );
            d(inin) = obj.in * d(inin);
            outout = ~inin & ~on;
            d(outout) = -obj.in * d(outout);
            % for nodes that were zero and are now non-zero, get the
            % sign by averaging the nodes of the touching elements that do
            % not contain the boundary
            ind = find( on & (abs(d)>obj.gerr) );
            for i1 = 1:length(ind)
                Ti = vertexAttachments(obj.mesh,ind(i1));
                Ni = obj.mesh.ConnectivityList( Ti{:}, : );
                Ni = unique( Ni(~any(ismember(Ni,C),2),:) );
                sig = sign(mean(obj.dist.Values(Ni(:))));
                d(ind(i1)) = sig * d(ind(i1));
            end
            obj = levelSet( X,'init', d, C );
        end
        % merge two meshes
        function dt = mergeMeshes(obj1,obj2)
            X = [ obj1.mesh.Points; obj2.mesh.Points ];
            N1 = size(obj1.mesh.Points,1);
            C = [ obj1.mesh.Constraints; obj2.mesh.Constraints+N1 ];
            % merge close-by nodes
            X = round( X/obj1.gerr ) * obj1.gerr;
            [X,~,ind] = unique( X, 'rows' );
            C = ind(C);
            if ~isempty(C)
                dt = delaunayTriangulation( X, C );
                dt = delaunayTriangulation( dt.Points, dt.Constraints );
            else
                dt = delaunayTriangulation( X );
            end
        end
        % reduce the constraints in one domain to be intersected
        function ind = reduceBoundary( obj, obj1 )
            out = ~inside( obj1, obj.mesh.Points, true );
            ind = ~all( out( obj.mesh.Constraints ), 2 );
        end
        % projection of a 2D domain onto a line
%         function line = projectLine( obj, x0, dir )
%             
%         end
     end
end
