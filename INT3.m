classdef INT3
    % definition of a mesh with INT3 elements
    %
    %  syntax: mesh = INT3( T, X )
    %
    %  X: [Nn*d] nodal coordinates matrix, where Nn is the number of nodes and
    %     d is the space dimension
    %  T: [Nel*3] element connectivity matrix indexing lines in X, where Nel is
    %     the number of elements and 3 the number of nodes in each element.
    %
    %  NB: the numeration of int3 elements is
    %           1---3---2
    %
    %  INT3 methods:
    %     plot               - Plot the TRI6 mesh with nodes and vertices
    %     elementsInBoundary - Returns the list of simplices inside a polygon
    %     subSet             - Extract a subset of the TRI6 mesh
    %     incenters          - middle of the segment
    %     freeBoundary       - Returns the nodes referenced by only one elemt
    %     parts              - separate the mesh into disjoint parts
    %     domain             - Returns a vector of size of the elements
    %     parts              - domain covered by the elements
    %     clean              - Removes unused X and repeated X and T
    %     size               - Number of elements
    %
    %  INT3 properties:
    %     T      - connectivity matrix [Ne*6]
    %     X      - nodal coordinates matrix [Nn*d]
    %     d      - space dimension
    %     Ne     - number of elements
    %     Nn     - number of nodes
    
    % R. Cottereau 07/2013
    
    properties
        Points           % Nn*1 coordinate matrix
        ConnectivityList % Ne*3 connectivity matrix
    end
    
    properties (Constant)
        gerr = 1e-9;  % error used to compute separations of points
    end
    
    properties (Dependent)
        Ne;     % number of intervals
        Nn;     % number of nodes
        Np;     % number of disjoint parts of the mesh
        d;      % space dimension
    end
    
    methods
        % define the INT3 object
        function obj = INT3( T, X, lclean )
            if nargin<3
                lclean = true;
            end
            obj.Points = X;
            obj.ConnectivityList = T;
            if lclean
                obj = clean( obj );
            end
        end
        % set the number of elements
        function Ne = get.Ne(obj)
            Ne = size(obj.ConnectivityList,1);
        end
        % set the number of nodes
        function Nn = get.Nn(obj)
            Nn = size(obj.Points,1);
        end
        % set the space dimension
        function d = get.d(obj)
            d = size(obj.Points,2);
        end
        % number of disjoint parts in the mesh
        function Np = get.Np(obj)
            Np = obj.Nn-obj.Ne;
        end
        % size
        function Ne = size(obj)
            Ne = size(obj.ConnectivityList);
        end
        % plot the mesh and the nodes
        function plot(obj)
            t = [obj.ConnectivityList(:,1);obj.ConnectivityList(end,2)];
            figure; plot( obj.Points(t), ones(size(t)), 'k-' );
            hold on; plot( obj.Points(t), ones(size(t)), 'ro' );
            bounds = max(obj.Points)-min(obj.Points);
            set(gca,'PlotBoxAspectRatio', [bounds 1 1] );
        end
        % Selects the elements of obj inside a given boundary
        function ind = elementsInBoundary( obj, ls, lall )
            if nargin<3
                lall = true;
            end
            ind = inside( ls, obj.Points, true );
            if lall
                ind = all( ind(obj.ConnectivityList), 2 );
            else
                ind = any( ind(obj.ConnectivityList), 2 );
            end
        end
        % returns a INT3 object using only a selected list of elements
        function [obj,ind] = subSet( obj, indT )
            obj = INT3( obj.ConnectivityList(indT,:), obj.Points, false );
            [obj,ind] = clean(obj);
        end
        % physical size of the elements
        function l = area(obj)
            l = sqrt( obj.Points(obj.ConnectivityList(:,2),:).^2 ...
                    - obj.Points(obj.ConnectivityList(:,1),:).^2 );
        end
        % get rid of nodes that are not used in T, and of repeated elements
        function [obj,indX] = cleanT(obj)
            % get rid of elements that are repeated
            obj.ConnectivityList = unique( obj.ConnectivityList, 'rows' );
            % get rid of elements that have area zero
            indT = area(obj)>0;
            obj.ConnectivityList = obj.ConnectivityList(indT,:);
            % get rid of nodes that are not used in T
            indX = unique( obj.ConnectivityList );
            for i1 = 1:length(indX)
                obj.ConnectivityList( obj.ConnectivityList==indX(i1) ) = i1;
            end
            obj.Points = obj.Points(indX,:);
        end
        % get rid of nodes that are repeated
        function [obj,indX] = cleanX(obj)
            xrnd = round(obj.Points/obj.gerr)*obj.gerr;
            [obj.Points,indX,indu] = unique( xrnd, 'rows' ,'first', 'legacy' );
            obj.ConnectivityList = reshape( indu(obj.ConnectivityList), obj.Ne, 2 );
        end
        % clean X of unused nodes and repeated nodes and elements
        function [obj,indX] = clean(obj)
            [obj,indX1] = cleanT(obj);
            [obj,indX2] = cleanX(obj);
            indX = indX1(indX2);
        end
        % centers of the segments
        function X = incenter(obj)
            X = (obj.Points(obj.ConnectivityList(:,1),:) ...
               + obj.Points(obj.ConnectivityList(:,2),:))/2;
        end
        % separate the mesh into disjoint parts
        function mi = parts(obj)
            mi = cell(obj.Np,1);
            dT = obj.ConnectivityList(2:end,1)-obj.ConnectivityList(1:end-1,2);
            dT = [1; find(dT~=0)+1; obj.Ne+1];
            for i1 = 1:obj.Np
                mi{i1} = subSet( obj, dT(i1):(dT(i1+1)-1) );
            end
        end
        % freeBoundary: nodes on the free boundary
        function [t,x,ind] = freeBoundary(obj)
            ind = sort(obj.ConnectivityList(:));
            ind = setdiff( ind, ind(diff(ind)==0) );
            t = (1:length(ind))';
            x = obj.Points(ind,:);
        end
        % levelSet describing the domain covered by the elements
        function ls = domain(obj)
            if obj.d>1
                error('not implemented yet')
            end
            Xb = zeros(obj.Np,2);
            mi = parts(obj);
            for i1 = 1:obj.Np
                Xb(i1,:) = [min(mi{i1}.Points) max(mi{i1}.Points)];
            end
            ls = levelSet1D( Xb(:,1), Xb(:,2), true );
        end
        % merge two meshes (find a mesh embedded in both obj1 and obj2)
        % it is not checked that both meshes are defined on same domain
        function obj = mergeMeshes( obj1, obj2, ~ )
            x = unique( [ obj1.Points; obj2.Points ] );
            N = length(x);
            t = [ (1:N-1)' (2:N)'];
            obj = INT3( t, x );
        end
        % bound an INT3 mesh inside a domain. Point are inserted at bnd
        function obj = bounded( obj, dom )
            [inX,on] = inside( dom, obj.Points, true );
            inT = all( inX(obj.ConnectivityList), 2 );
            bnd = any( inX(obj.ConnectivityList), 2 ) & ~inT ...
                & all( ~on(obj.ConnectivityList), 2 );
            if any(bnd)
                bnd1 = inX(obj.ConnectivityList(bnd,1));
                obj.ConnectivityList(bnd(bnd1),1) = obj.Nn+find(bnd1);
                obj.ConnectivityList(bnd(~bnd1),2) = obj.Nn+find(~bnd1);
                x1 = obj.Points(obj.ConnectivityList(bnd,1));
                x2 = obj.Points(obj.ConnectivityList(bnd,2));
                obj.Points = [ obj.Points; boundary( dom, x1, x2 ) ];
            end
            obj.ConnectivityList = obj.ConnectivityList(inT|bnd,:);
            obj = clean(obj);
        end
        % given the cartesian coordinates of points, return the barycentric
        function Xb = cartesianToBarycentric( obj, Tc, Xc )
            X1 = obj.Points(obj.ConnectivityList(Tc,1));
            X2 = obj.Points(obj.ConnectivityList(Tc,2));
            Xb = (Xc-X1)./(X2-X1);
        end
        % find elements in which nodes are located
        function ind = inElements( obj, Xc )
            ind = zeros(size(Xc,1),1);
            for i1 = 1:obj.Ne
                ind(Xc>=obj.Points(obj.ConnectivityList(i1,1)) ...
                   &Xc<=obj.Points(obj.ConnectivityList(i1,2))) = i1;
            end
        end
        % for each node in obj1, find the nodes that are inside the elts
        % that touch it, and the value of the linear FE basis function 
        % centered on Xr (=local barycentric coordinate of that node in the
        % element). Return a matrix in sparse format
        function [ Mx, My, Mval ] = XR2XI( obj, obj1 )
            ind = inElements( obj, obj1.Points );
            f1 = cartesianToBarycentric( obj, ind, obj1.Points );
            Mx = obj.ConnectivityList(ind,:);
            My = repmat((1:length(ind))',1,2);
            Mval = [1-f1 f1];
            ind = Mval>=obj.gerr;
            Mx = Mx(ind(:));
            My = My(ind(:));
            Mval = Mval(ind(:));
        end
    end
end