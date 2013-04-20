classdef TRI6
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

% R. Cottereau 03/2013

    properties
        X       % Nn*2 coordinate matrix
        T       % Ne*6 connectivity matrix
    end
    
    properties (Constant)
        gerr = 1e-9;  % error used to compute separations of points
    end
    
    properties (Dependent)
        d;      % space dimension
        Ne;     % number of elements
        Nn;     % number of nodes
        T3;     % underlying TRI3 elements
        X3;     % nodes of the underlying TRI3 elements
        ind3v6; % indices of TRI3 nodes in TRI6 numeration
        tri3;   % underlying mesh with TRI3 elements
    end

    methods
        % define the TRI6 object
        function obj = TRI6( T, X, lclean )
            if nargin<3
                lclean = true;
            end
            if size(T,2)==3
                ne = size(T,1);
                tri3 = TriRep(T,X);
                E = edges(tri3);
                x = tri3.X(:,1); y = tri3.X(:,2);
                X = [X; mean(x(E),2) mean(y(E),2)];
                T = [T T(:,1) zeros(ne,2)];
                for i1 = 3:-1:1
                    Xc = [x(T(:,i1))+x(T(:,i1+1)) y(T(:,i1))+y(T(:,i1+1))]/2;
                    [~,ind] = ismember( Xc, X, 'rows' );
                    T(:,i1+3) = ind;
                end
                obj = TRI6(T,X,lclean);
            elseif size(T,2)==6
                obj.X = X;
                obj.T = T;
                if lclean
                    obj = clean(obj);
                end
            else
                error('incorrect size of connectivity matrix T')
            end
        end
        % set the indices list of nodes in the TRI6 representation
        function ind = get.ind3v6(obj)
            ind = unique( obj.T(:,1:3), 'legacy' );
        end
        % set the element list of the underlying TriRep representation
        function T = get.T3(obj)
            T = obj.T(:,1:3);
            n3 = length(obj.ind3v6);
            for i1 = 1:n3
                T( T==obj.ind3v6(i1) ) = i1;
            end
        end
        % set the nodes of the TriRep representation
        function X = get.X3(obj)
            X = obj.X(obj.ind3v6,:);
        end
        % set the TriRep representation
        function tri3 = get.tri3(obj)
            tri3 = TriRep(obj.T3,obj.X3);
        end
        % set the number of elements
        function Ne = get.Ne(obj)
            Ne = size(obj.T,1);
        end
        % set the number of nodes
        function Nn = get.Nn(obj)
            Nn = size(obj.X,1);
        end
        % set the space dimension
        function d = get.d(obj)
            d = size(obj.X,2);
        end
        % returns a list of nodes that are not vertices
        function n = nodes(obj)
            n = unique( obj.T(:,4:6), 'legacy' );
        end
        % returns a list of vertices of the mesh
        function v = vertex(obj)
            v = unique( obj.T(:,1:3), 'legacy' );
        end
        % plot the mesh and the nodes
        function plot(obj)
            figure; triplot(obj.tri3,'color','k')
            hold on; scatter(obj.X(:,1),obj.X(:,2),50,'r','full');
        end
        % Selects the elements of obj inside a given boundary
        function ind = elementsInBoundary(obj,ls,lall)
            if nargin<3
                lall = true;
            end
            ind = inside( ls, obj.X, true );
            if lall
                ind = all( ind(obj.T), 2 );
            else
                ind = any( ind(obj.T), 2 );
            end
        end
        % returns a TRI6 object using only a selected list of elements
        function [obj,indX] = subSet(obj,indT)
            obj = TRI6( obj.T(indT,:), obj.X, false );
            [obj,indX] = clean(obj);
        end
        % get rid of nodes that are not used in T, and of repeated elements
        function [obj,indX] = cleanT(obj)
            % get rid of elements that are repeated
            [~,ind] = unique( sort(obj.T,2) ,'rows', 'first', 'legacy' );
            obj.T = obj.T(ind,:);
            % get rid of elements that have area zero
            xi = obj.X(:,1); yi = obj.X(:,2);
            ind = polyarea( xi(obj.T)', yi(obj.T)')'>obj.gerr;
            obj.T = obj.T(ind,:);
            % get rid of nodes that are not used in T
            indX = unique( obj.T, 'legacy' );
            for i1 = 1:length(indX)
                obj.T( obj.T==indX(i1) ) = i1;
            end
            obj.X = obj.X(indX,:);
        end
        % get rid of nodes that are repeated
        function [obj,indX] = cleanX(obj)
            xrnd = round(obj.X/obj.gerr)*obj.gerr;
            [obj.X,indX,indu] = unique( xrnd, 'rows' ,'first', 'legacy' );
            obj.T = indu(obj.T);
        end
        % clean X of unused nodes and repeated nodes and elements
        function [obj,indX] = clean(obj)
            [obj,indX1] = cleanT(obj);
            [obj,indX2] = cleanX(obj);
            indX = indX1(indX2);
        end
        % inherited from TriRep/size
        function N = size(obj)
            N = size(obj.tri3);
        end
        % inherited from TriRep/incenters
        function X = incenters(obj)
            X = incenters(obj.tri3);
        end
        % inherited from TriRep/freeBoundary
        function ls = freeBoundary(obj)
            [bnd,xf] = freeBoundary(obj.tri3);
            ls = levelSet( bnd, xf );
%             [~,xf] = freeBoundary(obj.tri3);
%             xf = [xf;xf(1,:)];
%             L = mean( max(obj.X3)-min(obj.X3) );
%             xm = min(obj.X3)-L/2;
%             xM = max(obj.X3)+L/2;
%             Xsup = [xm; [xm(1) xM(2)]; xM; [xM(1) xm(2)]];
%             dt = delaunayTriangulation( [obj.X3;Xsup] );
%             dt = TRI6( dt.ConnectivityList, dt.Points );
%             ls = levelSet( 'interface', xf, dt.X );
        end
    end
    
end