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
            if nargin<3
                lall = true;
            end
            if isa( ls, 'levelSet' )
                ind = inside( ls, obj.Points, true );
            elseif isa( ls, 'levelSet1D' )
                [proj,one2two] = projectLine( obj, ls.x0, ls.dir );
                ind = inside( ls, proj.Points, true );
                ind = ind(one2two(1:obj.Nn));
            end
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
            elseif isa( obj1, 'INT3' )
                % NB: one should implement 2nd order elements to be able to
                % compute these projections properly ... this
                % implementation is probably rather wrong !!!
                % computation of gradients in the y direction for all
                % elements of the TRI6
                [~,~,xloc,yloc] = projectLine( obj, obj1.x0, obj1.dir );
                xCut = obj1.Points;
                E = edges(obj);
                N = size(E,1);
                xc = cell(N,1);
                yc = cell(N,1);
                % list of cuts in the edges
                for i1 = 1:N
                    x1 = xloc(E(i1,1)); y1 = yloc(E(i1,1));
                    x2 = xloc(E(i1,2)); y2 = yloc(E(i1,2));
                    ind = find( xCut>=min(x1,x2) & xCut<=max(x1,x2) );
                    xc{i1} = ind;
                    yc{i1} = y1+(xCut(ind)-x1)/(x2-x1)*(y2-y1);
                    if isnan(yc{i1})
                        xc{i1} = [ind; ind];
                        yc{i1} = [y1; y2];
                    end
                end
                Mx = cell(obj.Ne,1);
                My = cell(obj.Ne,1);
                Mval = cell(obj.Ne,1);
                for i1 = 1:obj.Ne
                    Te = obj.ConnectivityList(i1,:);
                    Ei = sort([Te(1:2); Te(2:3); Te([3 1])],2);
                    lE = ismember( E, Ei, 'rows' );
                    xind = cat(1,xc{lE});
                    % exclude repeated cuts and cuts with one node only
                    xi = [xCut(xind) cat(1,yc{lE})];
                    xi = round( xi/obj1.gerr )*obj1.gerr;
                    [xi,ind] = unique( xi,'rows');
                    dx = diff(xi(:,1))==0;
                    dx = [false;dx] | [dx;false];
                    xi = xi(dx,:);
                    ind = ind(dx);
                    Nc = size(xi,1);
                    % position of the cut in local coordinates
                    bar = cartesianToBarycentric( obj, ...
                        i1*ones(Nc,1), xi )';
                    % length of the cut
                    h = diff(reshape(xi(:,2),2,Nc/2));
                    % should not count twice edges between elements
                    indj = any( bar==1, 1 );
                    indj = find(indj(:,1:2:end)&indj(:,2:2:end));
                    if size(indj,1)==0; indj = zeros(1,0); end
                    [indi,~] = find( bar(:,[2*indj-1;2*indj])==1 );
                    Ei = reshape(Te(indi),2,length(indj))';
                    nT = edgeAttachments( obj, Ei(:,1), Ei(:,2) );
                    nT = cellfun( @length, nT );
                    h(indj) = h(indj)./nT';
                    h = repmat( abs(h), [6 1]);
                    % y-derivative along the cut (assuming linear elements)
                    dy = [obj.Points(Te,:) ones(3,1)]\eye(3);
                    dy = repmat( dy(2,:)', [Nc/2 1] );
                    x1 = repmat( Te', [Nc 1] );
                    y1 = reshape(repmat(xind(ind),[1 3])',3*Nc,1);
                    val1 = h(:) .* bar(:)/2;
                    val2 = reshape(h(1:2:end,:),Nc/2*3,1) .* dy;
                    Mx{i1} = [x1; x1(1:end/2)];
                    My{i1} = [y1; y1(1:2:end)+obj1.Nn];
                    Mval{i1} = [val1; val2];
                end
                Mx = cat( 1, Mx{:} );
                My = cat( 1, My{:} );
                Mval = cat( 1, Mval{:} );
            else
                error(['XR2XI not implemented for this ' ...
                       'combination of dimensions'])
            end
            ind = abs(Mval(:))>obj.gerr;
            Mx = Mx(ind);
            Mval = Mval(ind);
            My = My(ind);
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
            E = edges(obj);
            % add intersection of edges that are cut by the line
            ind = yloc(E(:,1)).*yloc(E(:,2))<0;
            x1 = [xloc(E(ind,1)) yloc(E(ind,1))];
            x2 = [xloc(E(ind,2)) yloc(E(ind,2))];
            xc = x1(:,1)+x1(:,2).*(x2(:,1)-x1(:,1))./abs(x2(:,2)-x1(:,2));
            xloc = [ xloc; xc ];
            yloc = [ yloc; zeros(size(xc)) ];
            % sort and create INT3 object
            [xline,~,ind] = unique( xloc );
            N = length(xline);
            line = INT3( [(1:N-1)' (2:N)'], xline );
        end
    end
end