classdef TRI6 < TriRep
% definition of a mesh with TRI6 elements (extension of TriRep class to
% handle higher order representation of functions)
%
%  syntax: mesh = TRI6( T, X )
%
%  X: [Nn*d] nodal coordinates matrix, where d is the dimension of the
%     space considered, and Nn the number of nodes.
%  T: [Nel*t] element connectivity matrix indexing lines in X, where Nel is
%     the number of elements and t the number of nodes in each element.
%
%  Possible values of t are 3 and 6. If t=3, additional nodes are created
%  at the middle of all edges and and the corresponding mesh with TRI6
%  elements is created.
%
%  NB: the numeration of tri3 and tri6 elements is
%           1-------3       1---6---3
%            \     /         \     /
%             \   /           4   5
%              \ /             \ /
%               2               2
%
%  TRI6 methods:
%     plot               - Plot the TRI6 mesh with nodes and vertices
%     nodes              - Returns the list of nodes that are not vertices
%     vertex             - Returns the list of vertices
%     freeBoundary       - Returns the facets referenced by only one simplex
%     elementsInBoundary - Returns the list of simplices inside a polygon
%     subSet             - Extract a subset of the TRI6 mesh
%     cartToBary         - Converts the coordinates of a point from 
%                          cartesian to barycentric
%     size               - Returns the size of the Triangulation matrix
%     clean              - Removes unused X and repeated X and T
%
%  TRI6 inherited methods:
%     TRI6 does not inherits the methods of TriRep in a classical sense.
%     However, you can use all methods of TriRep with TRI6 and they will be
%     applied to the underlying TriRep class composed of TRI3 elements.
%     Refer to the help for TriRep for a list of these methods.
%
%  TRI6 properties:
%     T      - connectivity matrix [Ne*6]
%     X      - nodal coordinates matrix [Nn*d]
%     d      - space dimension
%     Ne     - number of elements
%     Nn     - number of nodes
%     tri3   - underlying mesh with TRI3 elements (TriRep class)
%     T3     - underlying TRI3 elements
%     X3     - nodes of the underlying tri3 mesh
%     ind3v6 - indices of the TRI3 nodes in the TRI6 nodes list

   % define the TRI6 object
    methods
        function obj = TRI6( T, X, lclean )
            obj = obj@TriRep(T,X(:,1),X(:,2));
            if nargin<3 || lclean
                obj = clean(obj);
            end
        end
    end
    
    properties (Constant)
        gerr = 1e-9;  % error used to compute separations of points
    end
    
    properties (Dependent)
        d;      % space dimension
        Ne;     % number of elements
        Nn;     % number of nodes
        DT;     % delaunay triangulation constrained
        DT2TRI; % index vector of the elements of DT into TRI6
    end

    methods
        % set the number of elements
        function Ne = get.Ne(obj)
            Ne = size(obj.Triangulation,1);
        end
        % set the number of nodes
        function Nn = get.Nn(obj)
            Nn = size(obj.X,1);
        end
        % set the space dimension
        function d = get.d(obj)
            d = size(obj.X,2);
        end
        % construct a constrained Delaunay triangulation attached to TRI6
        function DT = get.DT(obj)
            DT = DelaunayTri( obj.X, edges(obj) );
        end
        % indices of the DT elements in the TRI6 mesh
        function ind = get.DT2TRI(obj)
            T1 = sort( obj.Triangulation, 2 );
            [T1,ind1] = sortrows( T1, [1 2 3] );
            T2 = sort( obj.DT.Triangulation, 2 );
            [T2,ind2] = sortrows( T2, [1 2 3] );
            lind1 = ismember( T2, T1, 'rows' );
            lind2 = ismember( ind2, ind1, 'rows' );
            ind = NaN(size(ind2));
            ind(lind2) = ind1(ind2(ind2(lind1)));
        end
        % plot the mesh and the nodes
        function plot(obj)
            figure; triplot( obj, 'color', 'k' )
            hold on; scatter( obj.X(:,1), obj.X(:,2), 50, 'r', 'full' );
            Xc = incenters(obj);
            hold on; scatter( Xc(:,1), Xc(:,2), 50, 'bs' );
            bounds = max(obj.X) - min(obj.X);
            set( gca, 'PlotBoxAspectRatio', [bounds 1] );
        end
        % indices of elements containing specified points
        function TI = pointLocation(obj,qp)
            TI = pointLocation(obj.DT,qp);
            TI = obj.DT2TRI(TI);
        end
        % Selects the elements of obj inside a given boundary
        function ind = elementsInBoundary(obj,ls,lall)
            if nargin<3
                lall = true;
            end
            ind = inside( ls, obj.X, true );
            if lall
                ind = all( ind(obj.Triangulation), 2 );
            else
                ind = any( ind(obj.Triangulation), 2 );
            end
        end
        % returns a TRI6 object using only a selected list of elements
        function [obj,ind] = subSet(obj,indT)
            obj = TRI6( obj.Triangulation(indT,:), obj.X, false );
            [obj,ind] = clean(obj);
        end
        % domain covered by the mesh
        function ls = domain(obj)
            [bnd,xf] = freeBoundary(obj);
            ls = levelSet( bnd, xf );
        end
        % bound a mesh by a level-set. Nodes are added where the level-set
        % crosses elements
        function obj = bounded( obj, LSet )
            C = [obj.Triangulation(:,1:2); obj.Triangulation(:,2:3); ...
                 obj.Triangulation(:,[3 1])];
            Xm = obj.X;
            for i1=1:LSet.N
                C = [ C; LSet.T{i1}+size(Xm,1) ];
                Xm = [ Xm; LSet.X{i1} ];
            end
            [~,ind] = unique( sort(C,2) ,'rows', 'first');
            C = C(ind,:);
            ldt = clean( levelSet( C, Xm ) );
            obj = DelaunayTri( ldt.X{1}, ldt.T{1} );
            obj = TRI6( obj.Triangulation, obj.X );
            obj = subSet( obj, elementsInBoundary(obj,LSet,true) );
        end
        % merge two meshes (find a mesh embedded in both obj1 and obj2
        % and located within ls)
        function obj = MergeMeshes( obj1, obj2, ls )
            N = size(obj1.X,1);
            Xm = [ obj1.X ; obj2.X ];
            Tm = [ obj1.Triangulation ; obj2.Triangulation+N ];
            C = [Tm(:,1:2); Tm(:,2:3); Tm(:,[3 1])];
            [~,ind] = unique( sort(C,2) ,'rows', 'first');
            C = C(ind,:);
            ldt = clean( levelSet( C, Xm ) );
            obj = DelaunayTri( ldt.X{1}, ldt.T{1} );
            obj = TRI6( obj.Triangulation, obj.X );
            obj = subSet( obj, elementsInBoundary(obj,ls,true) );
        end
        % for each node in meshi, find nodes that are inside the elements
        % that touch it, and the value of the linear FE basis function 
        % centered on Xr (=local barycentric coordinate of that node in the
        % element). Return a matrix in sparse format
        function [ Mx, My, Mval ] = XR2XI( obj, obj1 )
            indx = pointLocation( obj, obj1.X );
            Mval = cartToBary( obj, indx, obj1.X );
            My = repmat( (1:size(Mval,1))', [3 1] );
            Mx = obj.Triangulation(indx,:);
            ind = abs(Mval(:))>obj.gerr;
            Mx = Mx(ind);
            Mval = Mval(ind);
            My = My(ind);
        end
        % get rid of nodes that are not used in T, and of repeated elements
        function [obj,indX] = cleanT(obj)
            % get rid of elements that are repeated
            T = obj.Triangulation;
            [~,ind] = unique( sort(T,2) , 'rows', 'first', 'legacy' );
            T = T(ind,:);
            % get rid of elements that have area zero
            xi = obj.X(:,1); yi = obj.X(:,2);
            ind = polyarea( xi(T)', yi(T)' )' > obj.gerr;
            T = T(ind,:);
            % get rid of nodes that are not used in T
            indX = unique( T, 'legacy' );
            for i1 = 1:length(indX)
                T( T==indX(i1) ) = i1;
            end
            obj = TRI6( T, obj.X(indX,:), false );
        end
        % get rid of nodes that are repeated
        function [obj,indX] = cleanX(obj)
            xrnd = round(obj.X/obj.gerr)*obj.gerr;
            [X,indX,indu] = unique( xrnd, 'rows' ,'first', 'legacy' );
            T = indu(obj.Triangulation);
            obj = TRI6( T, X, false );
        end
        % clean X of unused nodes and repeated nodes and elements
        function [obj,indX] = clean(obj)
            [obj,indX1] = cleanT(obj);
            [obj,indX2] = cleanX(obj);
            indX = indX1(indX2);
        end
    end
end