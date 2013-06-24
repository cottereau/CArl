function C = CouplingOperator( couple, Int, Rep, opt, Nmod)
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

if find(strcmp('HomeFE',Nmod))
[ x, y, C ] = CouplingOperatorHomeFE( couple.operator, Int.mesh.tri3, opt );
[ x, y, C ] = find( Rep.M * sparse( x, y, C ) * Int.M' );
        

C = struct( 'x', x, 'y', y, 'val', C );

if (strcmp( couple.mediator.type, 'stochastic' ))
    [ x, y, Cs ] = CouplingOperatorHomeFE( 'L2', Int.mesh.tri3, opt );
    Cs = Rep.M * sparse( x , y, Cs );
    [ xtheta, ~, Ctheta ] = find( sum( Cs, 2 ) );
    [ xpsi, ~, Cpsi ] = find( sum( Cs * Int.M', 1 )' );
    C.xBCpsi = xpsi;
    C.BCpsi = Cpsi;
    C.xtheta = xtheta;
    C.Ctheta = Ctheta;
end

elseif find(strcmp('FE2D',Nmod))
[ x, y, C ] = CouplingOperatorFE2D( couple.operator, Int.mesh.tri3, opt );
RepM = zeros(2*size(Rep.M));
RepM(1:2:end-1,1:2:end-1) = Rep.M;
RepM(2:2:end,2:2:end) = Rep.M;
%RepM(2:2:end,1:2:end-1) = Rep.M;
%RepM(1:2:end-1,2:2:end) = Rep.M;
IntM = zeros(2*size(Int.M));
IntM(1:2:end-1,1:2:end-1) = Int.M;
IntM(2:2:end,2:2:end) = Int.M;
%IntM(1:2:end-1,2:2:end) = Int.M;
%IntM(2:2:end,1:2:end-1) = Int.M;
[ x, y, C ] = find( RepM * sparse( x, y, C ) * IntM' );
C = struct( 'x', x, 'y', y, 'val', C );

if (strcmp( couple.mediator.type, 'stochastic' ))
    [ x, y, Cs ] = CouplingOperatorFE2D( 'L2', Int.mesh.tri3, opt );
    Cs = RepM * sparse( x , y, Cs );
    [ xtheta, ~, Ctheta ] = find( sum( Cs, 2 ) );
    [ xpsi, ~, Cpsi ] = find( sum( Cs * IntM', 1 )' );
    C.xBCpsi = xpsi;
    C.BCpsi = Cpsi;
    C.xtheta = xtheta;
    C.Ctheta = Ctheta;
end

elseif find((strcmp('TimobeamFE',Nmod)))
    %couplage Timo 2D
    
    [ x, y, C ] = CouplingOperatorTimo( couple.operator, Int.mesh.tri3, opt );
            
            M1 = Rep.Mbeam;
RepM = zeros(3*size(M1,1),3*size(M1,2));
RepM(1:size(M1,1),1:size(M1,2)) = M1;
RepM(size(M1,1)+1:2*size(M1,1),size(M1,2)+1:2*size(M1,2)) = M1;
%RepM(1:size(Rep.M,1),2:2:end) = Rep.M;
%RepM(size(Rep.M,1)+1:2*size(Rep.M,1),2:2:end) = Rep.M;
RepM(2*size(M1,1)+1:3*size(M1,1),2*size(M1,2)+1:3*size(M1,2)) = M1;
%RepM(2*size(Rep.M,1)+1:3*size(Rep.M,1),2:2:end) = Rep.M;

IntM = zeros(2*size(Int.M,1),2*size(Int.M,2));
IntM(1:2:end-1,1:2:end-1) = Int.M;
%IntM(size(Int.M,1)+1:2*size(Int.M,1),1:2:end-1) = Int.M;
IntM(2:2:end,2:2:end) = Int.M;
%IntM(size(Int.M,1)+1:2*size(Int.M,1),2:2:end) = Int.M;
%IntM(2*size(Int.M,1)+1:3*size(Int.M,1),1:2:end-1) = Int.M;
%IntM(2*size(Int.M,1)+1:3*size(Int.M,1),2:2:end) = Int.M;
[ x, y, C ] = find( RepM * sparse( x, y, C ) * IntM' );
C = struct( 'x', x, 'y', y, 'val', C );

if (strcmp( couple.mediator.type, 'stochastic' ))
    [ x, y, Cs ] = CouplingOperatorTimo( 'L2', Int.mesh.tri3, opt );
    Cs = RepM * sparse( x , y, Cs );
    [ xtheta, ~, Ctheta ] = find( sum( Cs, 2 ) );
    [ xpsi, ~, Cpsi ] = find( sum( Cs * IntM', 1 )' );
    C.xBCpsi = xpsi;
    C.BCpsi = Cpsi;
    C.xtheta = xtheta;
    C.Ctheta = Ctheta;
end
    
    
    
    
    
    
    
end