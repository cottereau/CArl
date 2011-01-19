function [ x, y, C ] = CouplingOperatorZacFE( operator, Int )
% COUPLINGOPERATORZACFE to construct the Arlequin coupling matrix 
% by calling a home-made FE code using the finite element toolboox
% of Rasmus Anthin
%
% syntax: [ x, y, C ] = CouplingOperatorHomeFE( model coupling, n )
%
%  model: structured array containing the information relative to the
%         model, in particular:
%       - 'type': describes the type of representation used for the
%                 geometry. This is used for all interpolation purposes,
%                 in particular for the definition of weight functions
%                 Implemented: {'FE' 'discrete'}
%       - 'code': code to be used to construct the stiffness matrices.
%                 Implemented: {'Comsol' 'HomeFE'}
%       - 'mesh': array dependent on 'type'.
%  coupling: cell of structured array describing the coupling options, in
%            particular
%       - 'mesh': indicates the subdomain of the mesh that corresponds to
%                 the coupling domain
%       - 'mediator': definition of the mediator space: 'default' or an
%                     integer indicating the optional choices (not always
%                     available)
%       - 'operator': coupling operator 'H1' 'L2'. The definition of these
%                     spaces may be dependent on the type of coupling
%  n: indicates what model is being considered [1 or 2]
%
%  output: the format is that of sparse matrices. The coupling matrix
%          is such that, schematically:
%               CouplingMatrix( x, y ) = C
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% C. Zaccardi 04/2010

% constants
Xi = Int.X;
Nni = size(Xi,1);

m = struct( 'mesh', Int, ...
    'load', 0, ...
    'property', struct('yg', ones(Nni,1),'dp',ones(Nni,1)), ...
    'BC', []);

% description of the model
X = m.mesh.X ;
T = m.mesh.T ;

% constants
Ne = size( T, 1 ) ;
[ Nn d ] = size( X ) ;

% property
E = m.property.yg ;
C = m.property.dp ;

if d==1
    % preparation for the use of genmat
    % assumption : the number of the node is in order
    lbc = find(T(1:Ne-1,2)~=T(2:Ne,1)) ;
    NX1end = T(lbc,2) ;
    NX2bg = T(lbc+1,1) ;
    % property has to be constant per element
    Eg = E ;
    Eg(2:end) = (1/2)*(E(1:end-1)+E(2:end)) ;
    % mesh separation
    X1 = X(1:NX1end,1) ;
    Eg1 = Eg(1:NX1end,1) ;
    C1 = C(1:NX1end,1) ;
    X2 = X(NX2bg:end,1) ;
    Eg2 = Eg(NX2bg:end,1) ;
    C2 = C(NX2bg:end,1) ;
    [Ki1,Fi1]=genmat1(X1',Eg1',C1',zeros(1,NX1end)) ;
    [Ki2,Fi2]=genmat1(X2',Eg2',C2',zeros(1,Nn-NX2bg+1)) ;
    Kze = zeros(NX1end,Nn-NX2bg+1) ;
    Ki = [Ki1, Kze ; Kze', Ki2] ;
    Fi = [Fi1, Kze ; Kze', Fi2] ;
elseif d==2
    [Ki,Fi]=genmat2(X,T,E,C,zeros(Nn,1)) ;
end
[ x, y, K ] = find( Ki );
[ z, k, F ] = find( Fi );

% choice of the coupling operator
switch operator

    % L2 coupling
    case 'L2'
        x = z;
        y = k;
        C = F;

    % H1 coupling
    case 'H1'
        x = [ x; z ];
        y = [ y; k ];
        C = [ K; F ];

    % unknown coupling operator
    otherwise
        error('unknown coupling operator')
end
