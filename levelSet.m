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
% considering a polygone defined by nodes Xp, and connectivity Tp.
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
    
    properties (Dependent)
        d     % dimension of the space
    end
    
    properties (Constant)
        gerr = 1e-9; % error used to compute separations of points
        in = -1;     % value of the distance function inside the domain
    end
    
    methods
        % Definition of level-set
        function obj = levelSet(varargin)
            if nargin==0
                obj.mesh = delaunayTriangulation;
                obj.dist = scatteredInterpolant;
            elseif nargin==1
                x0 = [0 0;1 0;1 1];
                v0 = obj.in*inf(3,1);
                obj.mesh = delaunayTriangulation( x0 );
                obj.dist = scatteredInterpolant( x0, v0, 'nearest' );
            elseif nargin==2||nargin==3
                x = varargin{1};
                v = varargin{2};
                obj.mesh = delaunayTriangulation( x );
                obj.dist = scatteredInterpolant( x, v, 'linear' );
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
                        x = [ xon; x; x1; xout ];
                        v = obj.in * [ zeros(nc,1); v; x2; -x2*ones(nc,1) ];
                        c = [(1:nc-1)' (2:nc)'];
                        obj.dist = scatteredInterpolant( x, v, 'linear' );
                        obj.mesh = delaunayTriangulation( x, c );
                        
                    % RECTANGLE OR SQUARE
                    case {'square','rectangle'}
                        xm = min(x2)/2;
                        de = 0.2*xm;
                        if numel(x2)==1
                            x2 = [x2 x2];
                        end
                        xv = x1(1) + [-de 0 x2(1) x2(1)+de];
                        yv = x1(2) + [-de 0 x2(2) x2(2)+de];
                        [xv,yv] = ndgrid(xv,yv);
                        xc = x1(1) + [xm x2(1)-xm]';
                        yc = x1(2) + [xm x2(2)-xm]';
                        x = [ [xc yc]; x];
                        xp = [xv(:) yv(:)];
                        tp = [6 7;7 11;11 10;10 6];
                        obj = levelSet( x, 'polygon', xp, tp );
                        
                    % POLYGON
                    case 'polygon'
                        if abs(x2(2:end,1)-x2(1:end-1,2))>obj.gerr
                            error('to be checked')
                        end
                        de = 0.1;
                        e1 = x1(x2(:,1),:);
                        e2 = x1(x2(:,2),:);
                        dir = e2 - e1;
                        dir = de * [dir(:,2) -dir(:,1)];
                        xsup = [e1+dir;e1-dir;e2+dir;e2-dir];
                        x = [ x1; x; xsup ];
                        v = -obj.in * dppoli( x, x2, x1 );
                        xv = x1(x2(:,1),1);
                        yv = x1(x2(:,1),2);
                        ind = inpolygon( x(:,1), x(:,2), xv, yv );
                        v(ind) = -v(ind);
                        obj.dist = scatteredInterpolant( x, v, 'linear' );
                        obj.mesh = delaunayTriangulation( x, x2 );

                    % RE-INITIALIZATION
                    case 'init'
                        obj.dist = scatteredInterpolant( x, x1, 'linear' );
                        obj.mesh = delaunayTriangulation( x, x2 );

                    otherwise
                        error('unknown type of geometrical domain')
                end
            end
            if (nargin==1||nargin==3||nargin==5) && ~varargin{end}
               obj = complement(obj);
            end
        end
        % get the dimension of the space
        function d = get.d(obj)
            d = size( obj.mesh.Points, 2 );
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
            V = obj.dist.Values;
            figure; trisurf(T,X(:,1),X(:,2),V);
            [ff,xf] = boundary(obj);
            x = [xf(ff(:,1),1) xf(ff(:,2),1)];
            y = [xf(ff(:,1),2) xf(ff(:,2),2)];
            hold on; plot( x', y', 'k-', 'LineWidth', 2 );
            shading flat; view(2); colorbar
        end
        % get the interface (for which distance=0)
        function varargout = boundary(obj)
            C = obj.mesh.Constraints;
            % attached elements are on different sides of the interface (or
            % on the boundary of the convex hull)
            on = true(size(C,1),1);
            bnd = ~all( ismember(C,convexHull(obj.mesh)), 2 );
            if ~isempty(bnd)
                TE = edgeAttachments( obj.mesh, C(bnd,1), C(bnd,2) );
                TE = cat( 1, TE{:} );
                [inXc,onXc] = inside( obj, incenter(obj.mesh,TE(:)), false );
                N = size(TE,1);
                inXc = reshape( inXc, N, 2 );
                onXc = reshape( onXc, N, 2 );
                on(bnd) = ~( all(inXc,2) | all(~(inXc|onXc),2) );
            end
            C = C(on,:);
            % both nodes of the constraint are on the boundary
            [~,on] = inside( obj, obj.dist.Points(C(:),:), true );
            C = C( all( reshape(on,size(C,1),2), 2 ), : );
            % a node inside the domain is not on two segments
%             [sc,ic] = sort(C(:));
%             dn = diff([ -Inf; sc; Inf ]);
%             dn = ~dn(2:end)==0 & ~dn(1:end-1)==0;
%             C = C( all( reshape( dn(ic), size(C) ), 2 ), : );
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
        % intersecting two domains composed of several level-sets
        function obj = intersection( obj1, obj2 )
            X = [ obj1.dist.Points; obj2.dist.Points ];
            N1 = size(obj1.dist.Points,1);
            C = [ obj1.mesh.Constraints; obj2.mesh.Constraints+N1 ];
            dt = delaunayTriangulation( X, C );
            V = max( [obj1.dist(dt.Points) obj2.dist(dt.Points)], [], 2 );
            C = dt.Constraints;
            C = C( all(abs(V(C))<obj1.gerr,2), : );
            obj = levelSet( dt.Points, 'init', V, C );
            obj = clean(obj);
        end
        % union of two level-sets
        function obj = union( obj1, obj2 )
            X = [ obj1.dist.Points; obj2.dist.Points ];
            N1 = size(obj1.dist.Points,1);
            C = [ obj1.mesh.Constraints; obj2.mesh.Constraints+N1 ];
            dt = delaunayTriangulation( X, C );
            V = min( [obj1.dist(dt.Points) obj2.dist(dt.Points)], [], 2 );
            C = dt.Constraints;
            C = C( all(abs(V(C))<obj1.gerr,2), : );
            obj = levelSet( dt.Points, 'init', V, C );  
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
            inon = inside(obj,obj.mesh.Points,true);
            outon = ~inside(obj,obj.mesh.Points,false);
            if all( inon )
                obj = levelSet(true);
            elseif all( outon )
                obj = levelSet(false);
            else
                C = boundary(obj);
                d = dppoli( obj.dist.Points, C, obj.dist.Points );
                inin = inon & ~outon;
                d(inin) = obj.in * d(inin);
                outout = outon & ~inon;
                d(outout) = -obj.in * d(outout);
                % for nodes that were zero and are now non-zero, get the
                % sign by averaging the touching nodes
                ind = find( (inon&outon) & (abs(d)>obj.gerr) );
                TI = vertexAttachments( obj.mesh, ind );
                for i1 = 1:length(ind)
                    ni = obj.mesh.ConnectivityList(TI{i1}',:);
                    d(ind(i1)) = sign(mean(d(ni(:)))) * d(ind(i1));
                end
                obj = levelSet( obj.dist.Points,'init', d, C );
            end
        end
     end
end
