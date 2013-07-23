classdef levelSet < delaunayTriangulation
% definition of a delaunayTriangulation class with additional methods
%
% LS = levelSet('polygon',Xp,Tp) creates a level-set object 
% considering a polygon defined by nodes Xp, and connectivity Tp.
%
% LS = levelSet(X) creates a level-set object, from the data of a set of
% nodes X and assuming that they, taken one after the other, create a
% bounded polygon
%
% LS = levelSet('circle',Xc,R) creates a level-set object considering a
% domain inside a circular interface centered on Xc, with radius R.
%
% LS = levelSet('square',Xc,L) creates a level-set object considering a
% domain inside a square interface with corner Xc and side length L.
% LS = levelSet(X,'rectangle',Xc,L) is the same with L a vector of two
% coordinates
%
% LS = levelSet(true) corresponds to an empty domain and
% LS = levelSet(false) corresponds to the full space.
% LS = levelSet is the same as levelSet(true)
%
%    levelSet properties:
%        in    - indicates whether we consider the inside or outside domain
%                [logical]
%
%    levelSet methods:
%     intersection - intersection of two domains
%     union        - union of two domains
%     complement   - complement of one domain within another
%     inside       - check whether points are inside the domain
    
    properties
        in  % indicates whether inside really means inside
    end
    
    properties (Constant)
        gerr = 1e-9; % error used to compute separations of points
    end
    
    methods
        % Definition of level-set
        function obj = levelSet(varargin)
            if nargin==0 
                X = [];
                C = [];
                lin = true;
            elseif nargin==1 && islogical(varargin{1})
                X = [];
                C = [];
                lin = varargin{1};
            elseif nargin==1
                X = varargin{1};
                C = [];
                lin = true;
            else
                geom = varargin{1};
                x1 = varargin{2};
                x2 = varargin{3};
                if nargin>3
                    lin = varargin{4};
                else
                    lin = true;
                end
                switch geom
                    % CIRCLE
                    case 'circle'
                        nc = 1e3;
                        t = linspace( 0, 2*pi, nc )';
                        X = [x1(1)+x2*cos(t) x1(2)+x2*sin(t)];
                        C = [(1:nc-1)' (2:nc)'];
                        
                    % RECTANGLE OR SQUARE
                    case {'square','rectangle'}
                        if numel(x2)==1
                            x2 = [x2 x2];
                        end
                        X = [x1(1)+[0 x2(1) x2(1) 0]' ...
                             x1(2)+[0 0 x2(2) x2(2)]'];
                        C = [1 2;2 3;3 4;4 1];
                        
                    % POLYGON
                    case 'polygon'
                        X = x;
                        C = x2;
 
                    otherwise
                        error('unknown type of geometrical domain')
                end
            end
            obj = delaunayTriangulation( X, C );
            obj.in = lin;
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
                obj.in = ~obj.in;
            else
                obj = intersection( obj, complement(obj1) );
            end
        end
        % intersecting two domains
        function obj = intersection(obj1,obj2)
            keyboard
            dt = mergeMeshes( obj1, obj2 );
            obj1 = updateDistance( obj1, dt.Points );
            obj2 = updateDistance( obj2, dt.Points );
            V = max( [obj1.dist(dt.Points) obj2.dist(dt.Points)], [], 2 );
            obj = levelSet( dt.Points, 'init', V, dt.Constraints );
        end
        % union of two domains
        function obj = union( obj1, obj2 )
            dt = mergeMeshes( obj1, obj2 );
            obj1 = updateDistance( obj1, dt.Points );
            obj2 = updateDistance( obj2, dt.Points );
            V = min( [obj1.dist(dt.Points) obj2.dist(dt.Points)], [], 2 );
            obj = levelSet( dt.Points, 'init', V, dt.Constraints );
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
            dt = delaunayTriangulation( X, C );
            dt = delaunayTriangulation( dt.Points, dt.Constraints );
        end
     end
end
