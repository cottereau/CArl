function c = DefineClassicalCoupling( c, m1 )
% DEFINECLASSICALCOUPLING to help in the definition of the coupling
% operator
%  syntax: cout = DefineClassicalCoupling( cin )
%  
%  cin : structured array containing the description of the coupling, with
%        the fields
%       -'type' : either 'zoom' or 'join'
%       -'models' : 1*2 vector naming the models to be coupled (the number
%                  k indicates model{k}). When 'zoom' is used, the first
%                  model of this vector should be the coarse one.
%       -'mediator' : structured array for the definition of the mediator 
%                  space. The fields are
%             --'type': 'deterministic' or 'stochastic' ('deterministic by
%                     default)
%             --'support': 1 or 2. indicates which model should be used as
%                     the support for the physical basis functions. Note
%                     that usually the coarser support should be used.
%       -'operator' : type of coupling operator used ('H1' or 'L2').
%       -'LevelSet1'; 'LevelSet2' : definition of level-sets bounding the
%             coupling zone, given as a mesh in (X,T) format, with the
%             appropriate dimension
%       -'epsilon' : residual value of the non-proeminent model
%
%  cout : same as cin, with the following additional fields
%       -'weight1': structured array defining the relevant information for
%                 model{coupling.models(1)} and containing the fields
%            --'ext': external geometrical limit of the coupling zone given 
%                 as a mesh in (X,T) format, with the appropriate dimension
%            --'int': internal geometrical limit of the coupling zone given
%                 in the same format as ext. The curve int should be 
%                 located inside the curve ext.
%              NB: in 1D, the level sets should include +/-Inf to give
%              sense to 'interior' and 'exterior' definitions
%            --'value': value of the weight function inside the coupling
%                       zone, given as a vector of coefficients of a
%                       polynomial (a scalar for a constant weight, a 2*1
%                       vector for a linear weight, etc ...)
%            --'extvalue': value of the weight outside the exterior curve
%            --'intvalue': value of the weight inside the interior curve
%       -'weight2': same as 'weight1' for the other model. When not
%                   defined, weight2=weight1, except for the values that
%                   are taken such that (value2+value1=1), and int and ext
%                   that are inversed.

% R. Cottereau 31/01/2012

% default case
if ~isfield( c, 'type' )
    return
end

% zoom or join
switch c.type
    case 'zoom'
        w1ext = c.epsilon;
        l = LSinLS( c.LevelSet1, c.LevelSet2 );
        
    case 'join'
        w1ext = 0;
        l = MDinLS( m1.mesh.X, c.LevelSet1 );
end
if l
    w1LSint = c.LevelSet2;
    w1LSext = c.LevelSet1;
else
    w1LSint = c.LevelSet1;
    w1LSext = c.LevelSet2;
end

% creation of weight function for first model
c.weight1 = struct( 'value', [1-c.epsilon c.epsilon], ...
                    'extvalue', w1ext, ...
                    'intvalue', 1, ...
                    'int', w1LSint, ...
                    'ext', w1LSext );

% creation of weight function for second model
c.weight2 = c.weight1;
c.weight2.int = c.weight1.ext;
c.weight2.ext = c.weight1.int;
c.weight2.intvalue = 1-c.weight1.extvalue;
c.weight2.extvalue = 1-c.weight1.intvalue;

function l = MDinLS( X, LS )
% MDINLS to indicate, in the case of join, whether points X are all inside 
% levelset LS or not
l = true;
switch size(LS.T,2)
    
    case 0 % 1D problem - 0D level set
        if any((X-eps)>max(LS.X))||any((X+eps)<min(LS.X))
            l = false;
        end
        
    case 1 % 2D problem - 1D level set
        xls = LS.X(LS.T(:,1),1);
        yls = LS.X(LS.T(:,1),2);
        if ~all(inpolygon( X(:,1), X(:,2), xls, yls ))
            l = false;
        end
end


