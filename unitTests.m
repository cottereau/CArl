function unitTests
% launches a series of unit Tests

% Tests for function XR2XI
X1 = [0 0;2 0;0 1]; T1 = [1 2 3];
X2 = [0 0;1 0;0 1]; T2 = [1 2 3];
[ Mx, My, Mval ] = XR2XI( TRI6(T1,X1), TRI6(T2,X2) );
l1 = norm(sparse(Mx,My,Mval)-[1 0 .5;0 1 0;0 0 .5])==0;
displayUnitTest( 'Test XR2XI 1', l1 );

X1 = [0 0;1 0;0 1]; T1 = [1 2 3];
X2 = [0 0;1 0;0 1]; T2 = [1 2 3];
[ Mx, My, Mval ] = XR2XI( TRI6(T1,X1), TRI6(T2,X2) );
l1 = norm(sparse(Mx,My,Mval)-eye(3))==0;
displayUnitTest( 'Test XR2XI 2', l1 );

X1 = [0 0;2 0;0 1;1 1/2]; T1 = [1 2 4;1 4 3];
X2 = [0 0;1 0;0 1;1 1/2]; T2 = [1 2 3;1 4 3];
[ Mx, My, Mval ] = XR2XI( TRI6(T1,X1), TRI6(T2,X2) );
l1 = norm(sparse(Mx,My,Mval)-[1 0 .5 0;0 1 0 0;0 0 0 1;0 0 .5 0])==0;
displayUnitTest( 'Test XR2XI 3', l1 );

% function display
function displayUnitTest( msg, l1 )
msg = [ ' ' msg ' ' repmat( '.', [1 30-length(msg)]) ' : ' ];
if l1
    msg = [msg 'OK'];
else
    msg = [msg 'FAILED'];
end
disp(msg);