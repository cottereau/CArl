function [ x, y, C ] = CouplingOperatorComsol( operator, Int )
% CouplingOperatorComsol to construct the Arlequin coupling matrix 
% by calling comsol
%
% For now, only linear element is implemented.
%
% syntax: [ x, y, C ] = CouplingOperatorComsol( operator, Int )
%  operator: coupling operator 'H1' 'L2'. The definition of these
%            spaces may be dependent on the type of coupling
%  Int: structured array containing the information relative to the mesh (X
%  and T)
%
%  output: the format is that of sparse matrices. The coupling matrix
%          is such that, schematically:
%               CouplingMatrix( x, y ) = C
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

% C. Zaccardi 07/2010

% message
disp('The couplingoperatorcomsol part is not optimized')

flclear fem

% COMSOL version
clear vrsn
vrsn.name = 'COMSOL 3.5';
vrsn.ext = 'a';
vrsn.major = 0;
vrsn.build = 603;
vrsn.rcs = '$Name:  $';
vrsn.date = '$Date: 2008/12/03 17:02:19 $';
fem.version = vrsn;

% mesh
X = Int.X' ;
T = Int.T' ;
dom = ones(1,size(T,2)) ;
el = cell(1,0);
el{1} = struct('type','tri','elem',T,'dom',dom) ;
m = femmesh(X,el);
fem.mesh = m ;

% Application mode coupling L2
clear appl
appl.mode.class = 'Helmholtz';
appl.sshape = 2;
appl.assignsuffix = '_hzeq';
clear prop
prop.elemdefault='Lag1';
appl.prop = prop;
clear equ
equ.c = 0 ;
equ.f = 0 ;
equ.ind = 1;
appl.equ = equ;
fem.appl{1} = appl;
fem.frame = {'ref'};
fem.border = 1;
clear units;
units.basesystem = 'SI';
fem.units = units;

% ODE Settings
clear ode
clear units;
units.basesystem = 'SI';
ode.units = units;
fem.ode=ode;

% Multiphysics
fem=multiphysics(fem);

% Extend mesh
fem.xmesh=meshextend(fem);
nodes = xmeshinfo(fem.xmesh,'out','nodes') ;

% Stifness matrix for the L2 coupling
[KL2] = assemble(fem,'Out',{'K'}) ;
KL2 = KL2(nodes.dofs,nodes.dofs) ;

flclear fem

% COMSOL version
clear vrsn
vrsn.name = 'COMSOL 3.5';
vrsn.ext = 'a';
vrsn.major = 0;
vrsn.build = 603;
vrsn.rcs = '$Name:  $';
vrsn.date = '$Date: 2008/12/03 17:02:19 $';
fem.version = vrsn;

% mesh
X = Int.X' ;
T = Int.T' ;
dom = ones(1,size(T,2)) ;
el = cell(1,0);
el{1} = struct('type','tri','elem',T,'dom',dom) ;
m = femmesh(X,el);
fem.mesh = m ;

% Application mode H1part
clear appl
appl.mode.class = 'Helmholtz';
appl.sshape = 2;
appl.assignsuffix = '_hzeq';
clear prop
prop.elemdefault='Lag1';
appl.prop = prop;
clear equ
equ.a = 0 ;
equ.f = 0 ;
equ.ind = 1;
appl.equ = equ;
fem.appl{1} = appl;
fem.frame = {'ref'};
fem.border = 1;
clear units;
units.basesystem = 'SI';
fem.units = units;

% ODE Settings
clear ode
clear units;
units.basesystem = 'SI';
ode.units = units;
fem.ode=ode;

% Multiphysics
fem=multiphysics(fem);

% Extend mesh
fem.xmesh=meshextend(fem);
nodes = xmeshinfo(fem.xmesh,'out','nodes') ;

% Stifness matrix for the H1 coupling
[KHpart] = assemble(fem,'Out',{'K'}) ;
KHpart = KHpart(nodes.dofs,nodes.dofs) ;

[ x, y, K ] = find( KHpart );
[ z, k, F ] = find( KL2 );

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

