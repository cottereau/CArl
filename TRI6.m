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
        function obj = TRI6(T,X)
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
                obj = TRI6(T,X);
            elseif size(T,2)==6
                obj.X = X;
                obj.T = T;
                obj = clean(obj);
            else
                error('incorrect size of connectivity matrix T')
            end
        end
        function T = get.T3(obj)
            T = obj.T(:,1:3);
            n3 = length(obj.ind3v6);
            for i1 = 1:n3
                T( T==obj.ind3v6(i1) ) = i1;
            end
        end
        function ind = get.ind3v6(obj)
            ind = unique( obj.T(:,1:3) );
        end
        function X = get.X3(obj)
            X = obj.X(obj.ind3v6,:);
        end
        function tri3 = get.tri3(obj)
            tri3 = TriRep(obj.T3,obj.X3);
        end
        function Ne = get.Ne(obj)
            Ne = size(obj.T,1);
        end
        function Nn = get.Nn(obj)
            Nn = size(obj.X,1);
        end
        function d = get.d(obj)
            d = size(obj.X,2);
        end
        function X = incenters(obj)
            X = incenters(obj.tri3);
        end
        function n = nodes(obj)
            n = unique(obj.T(:,4:6));
        end
        function v = vertex(obj)
            v = unique(obj.T(:,1:3));
        end
        function plot(obj)
            figure; triplot(obj.tri3,'color','k')
            hold on; scatter(obj.X(:,1),obj.X(:,2),50,'r','full');
        end
        function varargout = freeBoundary(obj)
            if nargout==1
                bnd = freeBoundary(obj.tri3);
                bnd = obj.ind3v6(bnd);
                varargout{1} = bnd;
            elseif nargout==2
                [varargout{1},varargout{2}] = freeBoundary(obj.tri3);
            end
        end
        function ind = elementsInBoundary(obj,bndT,bndX)
            xp = bndX(:,1);
            yp = bndX(:,2);
            ind = inpolygon( obj.X(:,1), obj.X(:,2), xp(bndT(:)), yp(bndT(:)) );
            ind = all( ind(obj.T), 2 );
        end
        function obj2 = subSet(obj,indT)
            obj2 = TRI6( obj.T(indT,:), obj.X );
        end
        function N = size(obj)
            N = size(obj.tri3);
        end
        % get rid of nodes that are not used in T, and of repeated elements
        function obj = cleanT(obj)
            ind = unique(obj.T);
            n3 = length(ind);
            for i1 = 1:n3
                obj.T( obj.T==ind(i1) ) = i1;
            end
            obj.X = obj.X(ind,:);
            [~,ind] = unique( sort(obj.T,2) ,'rows' );
            obj.T = obj.T(ind,:);
        end
        % get rid of nodes that are repeated
        function obj = cleanX(obj)
            xrnd = round(obj.X/obj.gerr)*obj.gerr;
            [~,indx,indu] = unique( xrnd, 'rows', 'stable' );
            Nx = size(obj.X,1);
            if length(indx)<Nx
                obj.X = obj.X(indx,:);
                obj.T = indu(obj.T);
            end
        end
        % clean X of unused nodes and repeated nodes and elements
        function obj = clean(obj)
            obj = cleanX(obj);
            obj = cleanT(obj);
        end
    end
    
end