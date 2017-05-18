classdef discontinuous1D
% Class for discontinous fonctions, defined as continuous functions over
% different 1D subdomains
%
% vd = discontinuous1D(ls,x,v) defines a continuous function, defined as a
% linear by parts function over each domain defined as a levelSet1D class.
% The nodes x defining the function do not have to be contained in the
% interior of the domain.
% 
% vd = discontinuous1D(ls1,x1,v1,ls2,x2,v2,...,lsn,xn,vn) defines a 
% discontinuous function where the v1,v2,...,vn are defined over subdomains
% S1, S2, ..., SN. Note that all functions must be defined with 3 inputs.
%
% NB: no verification is made on the overlapping of subdomains
%
%  discontinuous1D properties:
%     x   - Ns*1 cell of subdomains (levelSet1D class)
%     val - Ns*1 cell of structured array (x,f)
%     Ns  - number of subdomains
%
%  discontinuous1D methods:
%     interp    - Returns the values of the function at a series of points
%     plot      - plot the function over data set X
%     addRegion - define a new subdomain and corresponding function
%
% see levelSet

% R. Cottereau 03/2013

    properties
        x    % Ns*1 cell of subdomains
        val  % Ns*1 cell of structured array (x,f)
    end
    
    properties (Dependent)
        Ns;      % number of subdomains
    end
    
    methods
        % Definition of discontinuous function
        function obj = discontinuous1D(varargin)
            if nargin==0 || varargin{1}.N==0
                obj.x = {};
                obj.val = {};
            elseif nargin>=3
                obj.x{1} = varargin{1};
                if numel(varargin{3})==1
                    f = varargin{3}*ones(size(varargin{2},1),1);
                else
                    f = varargin{3};
                end
                obj.val{1} = struct( 'x', varargin{2}, 'f', f );
            else
                error('incorrect number of input parameters')
            end
            if nargin>3 && mod(nargin,3)==0
                for i1 = 2:nargin/3
                    xi = varargin{(i1-1)*3+1};
                    fi = varargin{(i1-1)*3+2};
                    vali = varargin{(i1-1)*3+3};
                    obj = addRegion( obj, discontinuous1D( xi, fi, vali ) );
                end
            end
        end
        % setting of Ns property
        function N = get.Ns(obj)
            N = length(obj.x);
        end
        % plotting discontinuous function
        function plot( obj, X )
            figure; plot( X, interp( obj, X ), '+-' );
        end
        % identify the subdomains where nodes are located
        function ind = inside( obj, X )
            ind = zeros( size(X,1), 1 );
            for i1 = 1:obj.Ns
                ind( inside(obj.x{i1}, X ) ) = i1;
            end            
        end
        % interpolation: getting values of f inside each subdomain
        function [fv,indg] = interp( obj, X, indg )
            fv = zeros( size(X,1), 1 );
            if nargin<3
                indg = zeros(size(X,1),1);
                for i1 = 1:obj.Ns
                    ind = inside( obj.x{i1}, X );
                    indg(ind) = i1;
                    fv(ind) = interp1( obj.val{i1}.x, obj.val{i1}.f, ...
                                             X(ind,:),'linear', 'extrap');
                end
            else
                for i1 = 1:obj.Ns
                    ind = indg==i1;
                    fv(ind) = interp1( obj.val{i1}.x, obj.val{i1}.f, ...
                                             X(ind,:),'linear', 'extrap');
                end
            end
        end
        % addregion: define a new subdomain and corresponding function
        function obj = addRegion(varargin)
            obj = varargin{1};
            if nargin==2
                obj1 = varargin{2};
                if obj1.Ns==0 && obj.Ns==0
                    obj = discontinuous;
                elseif obj.Ns==0
                    obj = obj1;
                elseif obj1.Ns>0
                    for i1 = 1:obj1.Ns
                        if obj1.x{i1}.N>0
                            obj.val{obj.Ns+1} = obj1.val{i1};
                            obj.x{obj.Ns+1} = obj1.x{i1};
                        end
                    end
                end
            elseif nargin==4
                obj1 = discontinuous1D(varargin{2},varargin{3},varargin{4});
                obj = addRegion( obj, obj1 );
            else
                error('incorrect number of arguments')
            end
        end
        % multiplication of two functions (where only one function is
        % defined, the output is by default empty)
        function obj = times( obj1, obj2 )
            if isfloat(obj2)
                obj = obj1;
                for i1 = 1:obj.Ns
                    obj.val{i1}.f = obj.val{i1}.f .* obj2;
                end
            elseif isa(obj2,'discontinuous1D')
                obj = discontinuous1D;
                for i1 = 1:obj1.Ns
                    for i2 = 1:obj2.Ns
                        x12 = intersection( obj1.x{i1}, obj2.x{i2} );
                        if x12.N>1
                            x1 = obj1.val{i1}.x;
                            x1 = x1( inside(obj1,x1)==i1, 1 );
                            x2 = obj2.val{i2}.x;
                            x2 = x2( inside(obj2,x2)==i2, 1 );
                            xt = [x1;x2];
                            f12 = f(obj1.val{i1},xt) .* f(obj2.val{i2},xt);
                            obj = addRegion( obj, x12, xt, f12 );
                        end
                    end
                end
            else
                error('discontinuous1D/times: not implemented yet')
            end
        end
        % power function
        function obj = power( obj1, n )
            obj = obj1;
            for i1=1:obj.Ns
                obj.val{i1}.f = (obj.val{i1}.f).^n;
            end
        end
        % addition of two functions
        function obj = plus( obj1, obj2 )
            if isfloat(obj2)
                obj = obj1;
                for i1 = 1:obj.Ns
                    obj.val{i1}.f = obj.val{i1}.f + obj2;
                end
            elseif isa(obj2,'discontinuous1D')
                obj = discontinuous1D;
                for i1 = 1:obj1.Ns
                    for i2 = 1:obj2.Ns
                        x12 = intersection( obj1.x{i1}, obj2.x{i2} );
                        if x12.N>1
                            x1 = obj1.val{i1}.x;
                            x1 = x1( inside(obj1,x1)==i1, 1 );
                            x2 = obj2.val{i2}.x;
                            x2 = x2( inside(obj2,x2)==i2, 1 );
                            xt = [x1;x2];
                            f12 = f(obj1.val{i1},xt) + f(obj2.val{i2},xt);
                            obj = addRegion( obj, x12, xt, f12 );
                        end
                    end
                end
            else
                error('discontinuous1D/plus: not implemented yet')
            end
        end
        % get full domain covered by the function
        function ls = fullDomain( obj )
            ls = levelSet1D( false );
            for i1 = 1:obj.Ns
                ls = union( ls, obj.x{i1} );
            end
        end
    end
end     
