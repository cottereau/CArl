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
        X       % Nn*1 coordinate matrix
        T       % Ne*3 connectivity matrix
    end
    
    properties (Constant)
        gerr = 1e-9;  % error used to compute separations of points
    end
    
    properties (Dependent)
        Ne;     % number of intervals
        Nn;     % number of nodes
        Nv;     % number of vertices
        Np;     % number of disjoint parts of the mesh
        d;      % space dimension
    end
    
    methods
        % define the INT3 object
        function obj = INT3( T, X, lclean )
            if nargin<3
                lclean = true;
            end
            if size(T,2)==2
                Xi = (X(T(:,1),:)+X(T(:,2),:))/2;
                T = [ T size(X,1)+(1:size(Xi,1))' ];
                obj = INT3( T, [X;Xi], lclean );
            elseif size(T,2)==3
                obj.X = X;
                obj.T = T;
                if lclean
                    obj = clean( obj );
                end
            else
                error('incorrect size of connectivity matrix T')
            end
        end
        % set the number of elements
        function Ne = get.Ne(obj)
            Ne = size(obj.T,1);
        end
        % set the number of nodes
        function Nn = get.Nn(obj)
            Nn = size(obj.X,1);
        end
        % set the number of vertices
        function Nv = get.Nv(obj)
            Nv = length(unique(obj.T(:,1:2)));
        end 
        % set the space dimension
        function d = get.d(obj)
            d = size(obj.X,2);
        end
        % number of disjoint parts in the mesh
        function Np = get.Np(obj)
            Np = obj.Nv-obj.Ne;
        end
        % size
        function Ne = size(obj)
            Ne = size(obj.T);
        end
        % plot the mesh and the nodes
        function plot(obj)
            figure;
            for i1 = 1:obj.Ne
                hold on;
                plot( obj.X(obj.T(i1,:),1), obj.X(obj.T(i1,:),2), 'k-' );
            end
            bounds = max(obj.X)-min(obj.X);
            set(gca,'PlotBoxAspectRatio', [bounds 1] );
        end
        % Selects the elements of obj inside a given boundary
        function ind = elementsInBoundary( obj, ls, lall )
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
        % returns a INT3 object using only a selected list of elements
        function obj = subSet( obj, indT )
            obj = INT3( obj.T(indT,:), obj.X );
        end
        % physical size of the elements
        function l = area(obj)
            l = sqrt( obj.X(obj.T(:,2),:).^2 - obj.X(obj.T(:,1),:).^2 );
        end
        % get rid of nodes that are not used in T, and of repeated elements
        function [obj,indX] = cleanT(obj)
            % get rid of elements that are repeated
            obj.T = unique( obj.T, 'rows' );
            % get rid of elements that have area zero
            indT = area(obj)>0;
            obj.T = obj.T(indT,:);
            % get rid of nodes that are not used in T
            indX = unique( obj.T );
            for i1 = 1:length(indX)
                obj.T( obj.T==indX(i1) ) = i1;
            end
            obj.X = obj.X(indX,:);
        end
        % get rid of nodes that are repeated
        function [obj,indX] = cleanX(obj)
            xrnd = round(obj.X/obj.gerr)*obj.gerr;
            [obj.X,indX,indu] = unique( xrnd, 'rows' ,'first', 'legacy' );
            obj.T = reshape( indu(obj.T), obj.Ne, 3 );
        end
        % clean X of unused nodes and repeated nodes and elements
        function [obj,indX] = clean(obj)
            [obj,indX1] = cleanT(obj);
            [obj,indX2] = cleanX(obj);
            indX = indX1(indX2);
        end
        % centers of the segments
        function X = incenters(obj)
            X = obj.X(obj.T(:,3),:);
        end
        % separate the mesh into disjoint parts
        function mi = parts(obj)
            Np = obj.Np;
            mi = cell(obj.Np,1);
            for i1 = 1:Np
                [~,~,bnd] = freeBoundary(obj);
                indT = bnd(1);
                tmp = [];
                while length(tmp)<length(indT)
                    tmp = indT;
                    indT = ismember( obj.T(:,1), obj.T(tmp,:) ) | ...
                           ismember( obj.T(:,2), obj.T(tmp,:) );
                end
                mi{i1} = subSet( obj, indT );
                obj = subSet( obj, ~indT );
            end
        end
        % freeBoundary: nodes on the free boundary
        function [t,x,ind] = freeBoundary(obj)
            t = obj.T(:,1:2);
            ind = sort(t(:));
            ind = setdiff( ind, ind(diff(ind)==0) );
            t = (1:length(ind))';
            x = obj.X(ind,:);
        end
        % levelSet describing the domain covered by the elements
        function ls = domain(obj)
            if obj.d>1
                error('not implemented yet')
            end
            Xb = zeros(obj.Np,1);
            mi = parts(obj);
            for i1 = 1:obj.Np
                Xb(i1,:) = [min(mi{i1}.X) max(mi{i1}.X)];
            end
            ls = levelSet1D( Xb(:,1), Xb(:,2), true );
        end
    end
end