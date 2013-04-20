classdef levelSet
% definition of a levelSet class for implicit definition of boundaries of
% complex convex or concave domains
%
% LS = levelSet(f) creates a level-set object, from the data of a 
% scatteredInterpolant function f, that is a distance function from an
% interface.
%
% LS = levelSet(x,fx) creates a level-set object, from the data of a set of
% nodes x and the value of the distance function fx evaluated at these 
% nodes.
%
% LS = levelSet('circle',Xc,R) creates a level-set object considering a
% domain inside a circular interface centered on Xc, with radius R.
%
% LS = levelSet('square',Xc,L) creates a level-set object considering a 
% domain inside a square interface with corner Xc and side length L.
%
% LS = levelSet('interface',Xi,Xv) creates a level-set object considering a 
% polygonal interface defined by nodes Xi, and evaluates the distance
% function at points Xv.
%
% LS = levelSet(...,in) considers the exterior of the convex hull of the 
% domain to be part of the domain if in==false, and the opposite if 
% in==true (default).
%
% LS = levelSet(true) corresponds to an empty domain and 
% LS = levelSet(false) corresponds to the full space.
%
%  levelSet properties:
%     dist - distance function from interface [scatteredInterpolant class]
%     in   - indicates how the function is extrapolated outside the convex
%            hull of dist [logical]
%
%  levelSet methods:
%     refine       - add new points to the definition of distance function
%     interface    - retrieve the interface as a cell of points lists
%     initialize   - recompute the distance function based on interface
%     intersection - intersection of two domains defined by levelsets
%     union        - union of two domains defined by levelsets
%     complement   - complement of a domain with respect to another
%     inside       - check whether points are in the domain
%     clean        - ???
%     checkConnex  - ???
%     plot         - plot the distance function and interface

% R. Cottereau 04/2013

    properties
        dist  % distance function from interface [scatteredInterpolant]
        in    % indicates how the function is extrapolated outside the convex
              % hull of dist [logical]
    end

    properties (Constant)
        gerr = 1e-9; % error used to compute separations of points
        discC = 1e3; % discretization for approximation of circle convexhull
        sign = -1;   % sign of distance function inside the interface
        Nx = 1e2;    % discretization for computation of contours
    end
    
    properties (Dependent)
        interface; % interface (with contourc format)
        Ni;        % number of disjoint interfaces
    end 
    
    methods
        % Definition of level-set
        function obj = levelSet(varargin)
            if nargin==0
                obj.dist = levelSet(true);
            elseif nargin==1 && islogical(varargin{1})
                obj.dist = scatteredInterpolant();
                obj.in = varargin{1};
            elseif isa(varargin{1},'scatteredInterpolant')
                obj.dist = varargin{1};
                if nargin>1
                    obj.in = varargin{2};
                else
                    obj.in = true;
                end
            % special case of ellipse/circle
            elseif strcmp(varargin{1},'circle')
                theta = linspace(0,2*pi,obj.discC)';
                R = varargin{3};
                Xc = varargin{2};
                X = [ Xc; [Xc(1)+R*cos(theta) Xc(2)+R*sin(theta)]];
                X = [ X; [Xc(1)+2*R*cos(theta) Xc(2)+2*R*sin(theta)]];
                F = obj.sign*[ R; zeros(obj.discC,1); -R*ones(obj.discC,1) ];
                obj.dist = scatteredInterpolant( X, F );
                if nargin>3
                    obj.in = varargin{4};
                else
                    obj.in = true;
                end
            % special case of rectangle/square
            elseif strcmp(varargin{1},'square')||strcmp(varargin{1},'rectangle')
                L = varargin{3};
                if length(L)==1
                    L = [L L];
                end
                % list of points
                Xc = varargin{2};
                de = 2*obj.gerr;
                de = 0.1;
                xv = Xc(:,1) + linspace(-.5,1.5,20)*L(1);
                xv = [xv Xc(:,1)+[-de de L(1)-de L(1)+de] ];
                yv = Xc(:,2) + linspace(-.5,1.5,20)*L(2);
                yv = [yv Xc(:,2)+[-de de L(2)-de L(2)+de] ];
                [xv,yv] = ndgrid(xv,yv);
                Xv = [xv(:) yv(:)];
                % four corners
                poli = [Xc;Xc+[L(1) 0];Xc+L;Xc+[0 L(2)];Xc];
                d = dppoli( Xv, poli );
                ind = inpolygon( Xv(:,1), Xv(:,2), poli(:,1), poli(:,2) );
                d(ind) = obj.sign*d(ind);
                d(~ind) = -obj.sign*d(~ind);
                obj.dist = scatteredInterpolant( Xv, d );
                if nargin>3
                    obj.in = varargin{4};
                else
                    obj.in = true;
                end
            % special case of interface definition
            elseif strcmp(varargin{1},'interface')
                Xv = varargin{3};
                Xp = varargin{2};
                Fv = dppoli( Xv, Xp );
                ind = inpolygon( Xv(:,1), Xv(:,2), Xp(:,1), Xp(:,2) );
                Fv(ind) = obj.sign*Fv(ind);
                Fv(~ind) = -obj.sign*Fv(~ind);
                obj = levelSet(Xv,Fv,true);
                if nargin>3
                    obj.in = varargin{4};
                else
                    obj.in = true;
                end
            else
                obj.dist = scatteredInterpolant(varargin{1},varargin{2});
                if nargin>2
                    obj.in = varargin{3};
                else
                    obj.in = true;
                end
            end
%            obj = clean(obj);
            obj.dist.ExtrapolationMethod = 'none';
        end
        % get the number of interface
        function Ni = get.Ni(obj)
            int = obj.interface;
            Ni = 0;
            n = 1;
            while n<size(int,1)
                n = n+int(n,2)+1;
                Ni = Ni+1;
            end
        end
        % get the definition of the interface
        function int = get.interface(obj)
            X = obj.dist.Points;
            xv = linspace( min(X(:,1)), max(X(:,1)), obj.Nx );
            yv = linspace( min(X(:,2)), max(X(:,2)), obj.Nx );
            [x,y] = ndgrid(xv,yv);
            V = reshape( obj.dist([x(:) y(:)]), obj.Nx, obj.Nx )';
            int = contourc( xv, yv, V, [0 0])';
%             % create mesh
%             dt = delaunayTriangulation( obj.dist.Points );
%             E = edges(dt);
%             V = obj.dist.Values;
%             X = obj.dist.Points;
%             % take all points on the interface inside
%             ind0 = abs(V)<obj.gerr;
%             V(ind0) = obj.sign*2*obj.gerr;
%             % get edges crossed by interface and touching elements
%             % get next possible nodes and edges crossed by level set
%             ind = find( V(E(:,1)).*V(E(:,2)) < -obj.gerr );
%             Ec = E(ind,:);
%             alpha = abs(V(Ec(:,1)))./abs(V(Ec(:,1))-V(Ec(:,2)));
%             alpha = repmat( alpha, [1 2] );
%             Xc = X(Ec(:,1),:) + alpha.*(X(Ec(:,2),:)-X(Ec(:,1),:));
%             Tn = edgeAttachments( dt, Ec );
%             nn = size(Tn,1);
%             Sn = zeros(nn,4);
%             for i1 = 1:nn
%                 Nn = setdiff(dt.ConnectivityList(Tn{i1},:),Ec(i1,:))';
%                 Sn(i1,:) = find(any(ismember(E,Nn),2) & any(ismember(E,Ec(i1,:)),2))';
%             end
%             Sn = Sn';
%             Sn = reshape( Sn( ismember(Sn(:),ind) ), 2, nn )' ;
%             Seg = struct( 'ind', find(ind), 'Xc', Xc, 'nextSeg', Sn );
%             % get single nodes that are on the interface
%             % get next possible nodes and edges crossed by level set
%             ind0 = find( abs(V)<obj.gerr );
%             nn = length(ind0);
%             Xc = X(ind0,:);
%             Tn = vertexAttachments( dt, ind0 );
%             Nn = inf(nn,10);
%             Sn = inf(nn,10);
%             for i1 = 1:nn
%                 Nni = setdiff(T(Tn{i1},:),ind0(i1));
%                 Nn(i1,1:length(Nni)) = Nni';
%                 Sni = find( all( ismember(E,Nni), 2 ));
%                 Sn(i1,1:length(Sni)) = Sni';
%             end
%             indall = T( all( abs(V(T))<obj.gerr, 2 ), : );
%             if ~isempty(indall)
%                 keyboard
%                 Nall = numel(indall);
%                 np = zeros(Nall,1);
%                 ns = zeros(Nall,1);
%                 for i1 = 1:Nall
%                     np(i1) = nnz(Nn(~ismember(ind0,indall),:)==indall(i1));
%                     ns(i1) = nnz(S.nextP(~ismember(ind0,indall),:)==indall(i1));
%                 end
%                 nindall = ~ismember(ind0,indall(n==0));
%                 ind0 = ind0( nindall );
%                 Nn = Nn( nindall, : );
%                 Xc = Xc( nindall, : );
%                 Sn = Sn( nindall, : );
%             end
%             P = struct( 'ind', ind0, 'Xc', Xc, 'nextP', Nn, 'nextSeg', Sn );
            % get rid of all-zero elements ???
            % get a series of single interface
%            int = cell(0,1);
%            while ~isempty(P.ind) || ~isempty(Seg.ind)
% while ~isempty(Seg.ind)
%                 [int,P,S] = singleInterface( obj, Seg );
% %                int{end+1,1} = p;
%             end
        end
        % get one interface only
        function ci = getInterface( obj, i1 )
            int = obj.interface;
            Ni = 0;
            n = 1;
            ci = [];
            while n<size(int,1)
                Ni = Ni+1;
                if Ni==i1
                    ci = int( n+(1:int(n,2)), : );
                    return
                end
                n = n+int(n,2)+1;
            end
        end
%         % get a single interface
%         function [int,P,S] = singleInterface( obj, S )
%             keyboard
%             ind = S.ind;
%             nextS = S.nextSeg(1,1);
%             for i1 = 2:size(ind,1)
%                 i2 = ind==nextS;
%                 tmp = ind( i2 );
%                 ind( i2 ) = ind( i1 );
%                 ind( i1 ) = tmp;
%                 nextS = setdiff( S.nextSeg(i2,:), ind(i1-1) )
%         %        i2 = find( S.nextSeg(ind(i1-1))
%             end
%            nextP = [];
%            nextS = [];
%            int = [];
%            if ~isempty(S.ind)
%                nextS = S.ind(1);
%            else
%                nextP = P.ind(1);
%            end
%             while ~isempty(P.ind) || ~isempty(S.ind)
%                 indS = find( ismember(S.ind,nextS), 1 );
%                 indP = find( ismember(P.ind,nextP), 1 );
%                 if ~isempty(indS)
%                     int = [int; S.Xc(indS,:)];
%                     nextP = S.nextP(indS,:);
%                     nextS = S.nextSeg(indS,:);
%                     nindS = [1:indS-1 indS+1:length(S.ind)];
%                     P.nextSeg( P.nextSeg==S.ind(indS) ) = inf;
%                     S.nextSeg( S.nextSeg==S.ind(indS) ) = inf;
%                     S.ind = S.ind(nindS);
%                     S.Xc = S.Xc(nindS,:);
%                     S.nextSeg = S.nextSeg(nindS,:);
%                     S.nextP = S.nextP(nindS,:);
%                 elseif ~isempty(indP)
%                     int = [int; P.Xc(indP,:)];
%                     nextP = P.nextP(indP,:);
%                     nextS = P.nextSeg(indP,:);
%                     nindP = [1:indP-1 indP+1:length(P.ind)];
%                     P.nextSeg( P.nextSeg==P.ind(indP) ) = inf;
%                     S.nextSeg( S.nextSeg==P.ind(indP) ) = inf;
%                     P.ind = P.ind(nindP);
%                     P.Xc = P.Xc(nindP,:);
%                     P.nextSeg = P.nextSeg(nindP,:);
%                     P.nextP = P.nextP(nindP,:);
%                 else
%                     keyboard
%                     error('to be checked ... does it only mean we closed the loop?')
%                 end           
%            end
%             int = [ int; int(1,:)];
%         end
        % re-compute the distance from the interface
        function obj = initialize( obj )
            keyboard
            X = obj.dist.Points;
            d = inf( size(X,1), 1 );
            for i1 = 1:obj.Ni
                Xp = getInterface( obj, i1 );
                Xp = [Xp;Xp(1,:)];
                di = dppoli( X, Xp );
                % change the sign for nodes inside the polygon
                ind = inpolygon( X(:,1), X(:,2), Xp(:,1), Xp(:,2));
                di(ind) = obj.sign*di(ind);
                di(~ind) = -obj.sign*di(~ind);
                [~,indj] = min(abs([d di]),[],2);
                d(indj==2) = di(indj==2);
            end
            obj = levelSet( X, d, obj.in );
        end
        % add points to the definition of a level-set
        function obj = refine(obj,X)
            X = unique( round(X/obj.gerr)*obj.gerr, 'rows' );
            V = obj.dist(X);
%            Xo = obj.dist.Points;
%            K = convhull(Xo);
%            ind = ~inpolygon( X(:,1), X(:,2), Xo(K,1), Xo(K,2) ) & (V==0);
            V(isnan(V)) = (-obj.in*obj.sign+~obj.in*obj.sign)*inf;
            obj = levelSet( X, V, obj.in );
        end
        % intersecting two domains composed of several level-sets
        function obj = intersection( obj1, obj2 )
            if ~isempty(obj1.dist)
                X1 = obj1.dist.Points;
            else
                X1 = [];
            end
            if ~isempty(obj2.dist)
                X2 = obj2.dist.Points;
            else
                X2 = [];
            end
            % obj1 is full space
            if isempty(X1)&&~obj1.in
                obj = obj2;
                return
            % obj2 is full space
            elseif isempty(X2)&&~obj2.in 
                obj = obj1;
                return
            % obj1 or obj2 are empty space
            elseif (isempty(X1)&&obj1.in) || (isempty(X2)&&obj2.in)
                obj = levelSet(true);
                return
            end
            % general case
            X = [ X1; X2];
            obj1 = refine( obj1, X );
            obj2 = refine( obj2, X );
            d = max( obj1.dist.Values, obj2.dist.Values );
            obj = levelSet( obj1.dist.Points, d, true );
%            obj = initialize( obj );
        end
        % intersecting two domains composed of several level-sets
        function obj = union( obj1, obj2 )
%             % obj1 is empty space
%             if (obj1.Ni==0&&obj1.in)
%                 obj = obj2;
%                 return
%             % obj2 is empty space
%             elseif (obj2.Ni==0&&obj2.in)
%                 obj = obj1;
%                 return
%             % obj1 or obj2 are full space
%             elseif (obj1.Ni==0&&~obj1.in) || (obj2.Ni==0&&~obj2.in)
%                 obj = levelSet(false);
%                 return
%             end
            % general case
            X = [ obj1.dist.Points; obj2.dist.Points];
            obj1 = refine( obj1, X );
            obj2 = refine( obj2, X );
            d = min( obj1.dist.Values, obj2.dist.Values );
            obj = levelSet( obj1.dist.Points, d, true );
            obj = initialize( obj );
        end
        % complement of obj2 in obj1
        function obj = complement( obj1, obj2 )
            if nargin==1
                d = obj1.dist;
                obj = levelSet( d.Points, -d.Values, ~obj1.in );
            else
                obj = intersection( obj1, complement(obj2) );
            end
        end
        % check whether points are inside the level-set
        function l = inside( obj, X, lon )
            if nargin<3 || lon
                l = obj.dist(X)<=obj.gerr;
            else
                l = obj.dist(X)<=-obj.gerr;
            end
        end
        % clean X of unused and repeated nodes
        function obj = clean(obj)
            % get rid of infinity values
            ind = ~isinf(obj.dist.Values);
            X = obj.dist.Points(ind,:);
            V = obj.dist.Values(ind);
            % get rid of repeated nodes
            [X,ind] = unique( round(X/obj.gerr)*obj.gerr, 'rows' );
            V = V(ind);
            % create level set
            obj.dist = scatteredInterpolant( X, V );
        end
        % plot function
        function  plot( obj )
            figure;
            X = obj.dist.Points;
            dt = delaunayTriangulation(X);
            trisurf(dt.ConnectivityList,X(:,1),X(:,2),obj.dist.Values);
            for i1 = 1:obj.Ni
                ci = getInterface( obj, i1 );
                hold on; plot( ci(:,1), ci(:,2), 'k-', 'linewidth', 2 );
            end
            view(2);
        end
    end
end
