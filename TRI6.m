classdef TRI6 < triangulation
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
%  TRI6 is a subclass of triangulation and hence inherits its properties
%  and methods
%
%  TRI6 properties:
%     d      - space dimension
%     Ne     - number of elements
%     Nn     - number of nodes
%     DT     - delaunay triangulation constrained
%     DT2TRI - index vector of the elements of DT into TRI6

   % define the TRI6 object
    methods
        function obj = TRI6( T, X, lclean )
            obj = obj@triangulation(T,X(:,1),X(:,2));
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
        % construct a constrained Delaunay triangulation attached to TRI6
        function DT = get.DT(obj)
            DT = delaunayTriangulation( obj.Points, obj.edges );
        end
        % indices of the DT elements in the TRI6 mesh
        function ind = get.DT2TRI(obj)
            T1 = sort( obj.ConnectivityList, 2 );
            [T1,ind1] = sortrows( T1, [1 2 3] );
            T2 = sort( obj.DT.ConnectivityList, 2 );
            [T2,ind2] = sortrows( T2, [1 2 3] );
            lind1 = ismember( T2, T1, 'rows' );
            ind = NaN(size(ind2));
            ind(ind2(lind1)) = ind1;
        end
        % plot the mesh and the nodes
        function plot(obj)
            figure; triplot( obj, 'color', 'k' )
            hold on; scatter( obj.Points(:,1), obj.Points(:,2), 50, 'r', 'full' );
            Xc = incenter(obj);
            hold on; scatter( Xc(:,1), Xc(:,2), 50, 'bs' );
            bounds = max(obj.Points) - min(obj.Points);
            set( gca, 'PlotBoxAspectRatio', [bounds 1] );
        end
        % indices of elements containing specified points
        function TI = pointLocation(obj,qp)
            TI = pointLocation(obj.DT,qp);
            TI = obj.DT2TRI(TI);
        end
        % Selects the elements of obj inside a given boundary
        function ind = elementsInBoundary(obj,ls,lall)
            if ~isa( ls, 'levelSet' )
                ind = true(obj.Ne,1);
                return
            end
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
        % returns a TRI6 object using only a selected list of elements
        function [obj,ind] = subSet(obj,indT)
            obj = TRI6( obj.ConnectivityList(indT,:), obj.Points, false );
            [obj,ind] = clean(obj);
        end
        % domain covered by the mesh
        function ls = domain(obj)
            [bnd,xf] = freeBoundary(obj);
            ls = levelSet( obj.Points, 'polygon', xf, bnd );
        end
        % bound a mesh by a level-set. Nodes are added where the level-set
        % crosses elements
        function obj = bounded( obj, LSet )
            X = [ obj.Points; LSet.mesh.Points ];
            C = [ edges(obj); LSet.mesh.Constraints+size(obj.Points,1)];
            obj = delaunayTriangulation( X, C );
            obj = delaunayTriangulation( obj.Points, obj.Constraints );
            obj = TRI6( obj.ConnectivityList, obj.Points );
            [inin,on] = inside( LSet, obj.Points, false );
            indT = any(inin(obj.ConnectivityList),2) | ...
                   all(on(obj.ConnectivityList),2);
            obj = subSet( obj, indT );
        end
        % merge two meshes (find a mesh embedded in both obj1 and obj2
        % and located within ls)
        function obj = mergeMeshes( obj1, obj2 )
            if isa( obj2, 'TRI6' )
                X = [ obj1.Points ; obj2.Points ];
                C = [ edges(obj1); edges(obj2)+size(obj1.Points,1) ];
                obj = delaunayTriangulation( X, C );
                obj = TRI6( obj.ConnectivityList, obj.Points );
            elseif isa( obj2, 'INT3' )
                obj = mergeMeshes( obj2, obj1 );
            else
                error(['mergeMeshes not implemented for this ' ...
                       'combination of dimensions'])
            end
        end
        % for each node in meshi, find nodes that are inside the elements
        % that touch it, and the value of the linear FE basis function 
        % centered on Xr (=local barycentric coordinate of that node in the
        % element). Return a matrix in sparse format
        function [ Mx, My, Mval ] = XR2XI( obj, obj1 )
            if isa( obj1, 'TRI6' )
                indx = pointLocation( obj, obj1.Points );
                Mval = cartesianToBarycentric( obj, indx, obj1.Points );
                My = repmat( (1:size(Mval,1))', [3 1] );
                Mx = obj.ConnectivityList(indx,:);
                ind = abs(Mval(:))>obj.gerr;
                Mx = Mx(ind);
                Mval = Mval(ind);
                My = My(ind);
            elseif isa( obj1, 'INT3' )
                % NB: one should implement 2nd order elements to be able to
                % compute these projections properly ... this
                % implementation is probably rather wrong !!!
                % computation of gradients in the y direction for all
                % elements of the TRI6
                gry = zeros(obj.Ne,3);
                for i1 = 1:obj.Ne
                    Xe = obj.Points(obj.ConnectivityList(i1,:),:);
                    gr = [Xe ones(3,1)]\eye(3);
                    gry(i1,:) = gr(2,:);
                end
                Mx = cell(obj.Ne,1);
                My = cell(obj.Ne,1);
                Mval = cell(obj.Ne,1);
                [~,~,xloc,yloc] = projectLine( obj, [0 0], [1 0] );
                edg = edges( obj );
                for i1 = 1:obj1.Nn
                    % finding edges intersected by a section at x=xloc
                    edgCut = xloc(edg(:,1))<=obj1.Points(i1) ...
                                          & xloc(edg(:,2))>=obj1.Points(i1);
                    edgCut = edg(edgCut,:);
                    X1 = [xloc(edgCut(:,1)) yloc(edgCut(:,1))];
                    X2 = [xloc(edgCut(:,2)) yloc(edgCut(:,2))];
                    phi = (obj1.Points(i1)-X1(:,1)) ./ (X2(:,1)-X1(:,1));
                    yCut = X1(:,2)+ phi.*(X2(:,2)-X1(:,2));
                    inan = isnan(yCut);
                    yCut(inan) = (X1(inan,2)+X2(inan,2))/2;
                    [yCut,ind] = sort(yCut);
                    edgCut = edgCut(ind,:);
                    h = diff(yCut);
                    h = [[0;h] [h;0]];
                    h(find(inan(ind))-1,2) = 0;
                    h(find(inan(ind))+1,1) = 0;
                    h = sum(h,2)/2;
                    % computing average of phi
                    phi = [1-phi(ind,:) phi(ind,:)];
                    phi(inan(ind),:) = 1;
                    m1 = [ phi(:,1).*h/2; phi(:,2).*h/2 ];
                    % computing average of grad phi
                    TI = edgeAttachments( obj, edgCut(:,1), edgCut(:,2));
                    gri = zeros( size(edgCut,1),2);
                    for i2 = 1:size(gri,1);
                        tt = obj.ConnectivityList(TI{i2},:);
                        gri1 = gry( TI{i2}(1), tt(1,:)==edgCut(i2,1) );
                        gri2 = gry( TI{i2}(1), tt(1,:)==edgCut(i2,2) );
                        if size(tt,1)>1
                        gri1 = (gri1+gry(TI{i2}(2),tt(2,:)==edgCut(i2,1)))/2;
                        gri2 = (gri2+gry(TI{i2}(2),tt(2,:)==edgCut(i2,2)))/2;
                        end
                        gri(i2,:) = [gri1 gri2];
                    end
                    m2 = [ gri(:,1).*h/2; gri(:,2).*h/2 ];
                    % storing
                    Mval{i1} = [ m1; m2(:) ];
                    Mx{i1} = repmat( edgCut(:), [2 1] );
                    on = ones(size(Mx{i1},1)/2,1);
                    My{i1} = [i1*on; (i1+obj1.Nn)*on];
                    ind = abs(Mval{i1}(:))>obj.gerr;
                    Mx{i1} = Mx{i1}(ind);
                    My{i1} = My{i1}(ind);
                    Mval{i1} = Mval{i1}(ind);
                end
                Mx = cat(1,Mx{:});
                My = cat(1,My{:});
                Mval = cat(1,Mval{:});
            else
                error(['XR2XI not implemented for this ' ...
                       'combination of dimensions'])
            end
        end
        % get rid of nodes that are not used in T, and of repeated elements
        function [obj,indX] = cleanT(obj)
            % get rid of elements that are repeated
            T = obj.ConnectivityList;
            [~,ind] = unique( sort(T,2) , 'rows', 'first', 'legacy' );
            T = T(ind,:);
            % get rid of elements that have area zero
            xi = obj.Points(:,1); yi = obj.Points(:,2);
            ind = polyarea( xi(T)', yi(T)' )' > obj.gerr;
            T = T(ind,:);
            % get rid of nodes that are not used in T
            indX = unique( T, 'legacy' );
            for i1 = 1:length(indX)
                T( T==indX(i1) ) = i1;
            end
            obj = TRI6( T, obj.Points(indX,:), false );
        end
        % get rid of nodes that are repeated
        function [obj,indX] = cleanX(obj)
            xrnd = round(obj.Points/obj.gerr)*obj.gerr;
            [X,indX,indu] = unique( xrnd, 'rows' ,'first', 'legacy' );
            T = indu(obj.ConnectivityList);
            obj = TRI6( T, X, false );
        end
        % clean X of unused nodes and repeated nodes and elements
        function [obj,indX] = clean(obj)
            [obj,indX1] = cleanT(obj);
            [obj,indX2] = cleanX(obj);
            indX = indX1(indX2);
        end
        % projection of a 2D mesh onto a 1D mesh
        function [line,ind,xloc,yloc] = projectLine( obj, x0, dir )
            x = obj.Points(:,1)-x0(:,1);
            y = obj.Points(:,2)-x0(:,2);
            xloc = [x y]*dir';
            yloc = [x y]*[-dir(2); dir(1)];
            [xline,~,ind] = unique( xloc );
            N = length(xline);
            line = INT3( [(1:N-1)' (2:N)'], xline );
        end
    end
end