function C = CouplingOperator( couple, Int, Rep, opt )
% COUPLINGOPERATOR to construct the coupling operator
%
% syntax: [ x, y, C ] = CouplingOperator( model, coupling, n )
%
%  couple: structured array containing the fields
%       - 'code': code that should be used to compute the matrix
%       - 'operator': coupling operator 'H1' 'L2'. The definition of these
%                     spaces may be dependent on the type of coupling
%  Int: structured array containing the passage matrix from the coupling
%       matrix in integration DOFs into coupling DOFs
%  Rep: structured array containing the passage matrix from the coupling
%       matrix in integration DOFs into model DOFs
%
%  C is a sparse matrix
%
% copyright: Laboratoire MSSMat, Ecole Centrale Paris - CNRS UMR 8579
% contact: regis.cottereau@ecp.fr

switch lower(couple.code)
    
    % ACOUSTIC COUPLING
    case {'homefe','montecarlohomefe'}
        C = CouplingOperatorHomeFE( couple.operator, Int.mesh, opt );
        C = Rep.M * C * Int.M';
        
        if strcmpi( couple.code, 'montecarlohomefe' )
            Cs = Rep.M * CouplingOperatorHomeFE( 'L2', Int.mesh, opt );
            N = size(C,2);
            C = [ C sum( Cs, 2 ) sparse(size(Cs,1),1);
                  sparse(N,N+1) sum( Cs * Int.M', 1 )'];
        end
        
    % ELASTIC COUPLING
    case 'fe2d'
        C = CouplingOperatorFE2D( couple.operator, Int.mesh, opt );
        C = Rep.M * C * Int.M';
      
    % TIMOSCHENKO BEAM COUPLING
    case 'beam'
        opt.h = couple.h;
        C = CouplingOperatorBeam( couple.operator, Int.mesh, opt );
        C = Rep.M * C * Int.M';
              
    % UNKNOWN COUPLING CASE
    otherwise
        error('unknown coupling case')
        
end