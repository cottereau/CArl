function C = CouplingOperator( couple, Int, Rep, opt )
% COUPLINGOPERATOR to construct the coupling operator
% 
% syntax: [ x, y, C ] = CouplingOperator( model, coupling, n )
%
%  model: structured array containing the information relative to the
%         model, in particular:
%       - 'type': describes the type of representation used for the
%                 geometry. This is used for all interpolation purposes,
%                 in particular for the definition of weight functions
%                 Implemented: {'FE' 'discrete'}
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

% R. Cottereau 04/2010


[ x, y, C ] = CouplingOperatorHomeFE( couple.operator, Int.mesh, opt );
[ x, y, C ] = find( Rep.M * sparse( x, y, C ) * Int.M' );
C = struct( 'x', x, 'y', y, 'val', C );

if strcmp( couple.mediator.type, 'stochastic' )
    [ x, y, Cs ] = CouplingOperatorHomeFE( 'L2', Int.mesh, opt );
    Cs = Rep.M * sparse( x , y, Cs );
    [ xtheta, ~, Ctheta ] = find( sum( Cs, 2 ) );
    [ xpsi, ~, Cpsi ] = find( sum( Cs * Int.M', 1 )' );
    C.xBCpsi = xpsi;
    C.BCpsi = Cpsi;
    C.xtheta = xtheta;
    C.Ctheta = Ctheta;
end




