function [ x, y, C ] = CouplingOperatorTimo( operator, mesh, opt )

coordinates = mesh.X(mesh.X(:,2)==0,1);
numberNodes = size(coordinates,1);
elementNodes=[(1:numberNodes-1)' (2:numberNodes)'];
[Nne,ne] = size(elementNodes);

% computation of MassMatrix
m = struct( 'mesh', mesh, ...
            'property', ones(Nne,ne), ...
            'load', zeros(Nne,ne,3), ...
            'BC', [] );
        
numberElements = size(elementNodes,1);
GDof=3*numberNodes;
thickness=max(coordinates)-min(coordinates);
I=thickness^3/12;
for ijk = 1:numberElements
    Comp{ijk} = [ I 0 0; 0 1 0;0 0 thickness];
end
[stiffness,~,mass]=formStiffnessMassTimoshenkoBeam(GDof,numberElements,elementNodes,numberNodes,coordinates,Comp,m.load,1,I,thickness);

[x,y] = find(mass);
indtemp = (find(mass));
C = mass(indtemp);
%maxx = size(mass,1);
% choice of the coupling operator
switch operator
    
    % L2 coupling
    case 'L2'            
        
    % H1 coupling
    case 'H1'
       % [stiffC1,stiffC2] = stifness_coupling_Timo( elements,coordinates);
        [z,k] = find(stiffness);
        indtemp= find(stiffness);
        K2=(stiffness(indtemp));
        x = [ z; x ];
        y = [ k; y ];
        C = [ opt.kappa*K2; C ];
        
             
    % unknown coupling operator
    otherwise
        error('unknown coupling operator')
end





function [stiffness,force,mass]=...
formStiffnessMassTimoshenkoBeam(GDof,numberElements,...
elementNodes,numberNodes,xx,C,P,rho,I,thickness)
% computation of stiffness matrix and force vector
% for Timoshenko beam element
stiffness=zeros(GDof);
mass=zeros(GDof);
force=zeros(GDof,1);
% stiffness matrix
gaussLocations=[0.577350269189626;-0.577350269189626];
gaussWeights=ones(2,1);
%gaussLocations=[0.];
%gaussWeights=[2.];
% bending contribution for stiffness matrix
for e=1:numberElements
indice=elementNodes(e,:);
elementDof=[ indice+numberNodes indice+2*numberNodes];
indiceMass=indice+2*numberNodes;
ndof=length(indice);
length_element=xx(indice(2))-xx(indice(1));
detJacobian=length_element/2;invJacobian=1/detJacobian;
Ce = C{e};
for q=1:size(gaussWeights,1) ;
pt=gaussLocations(q,:);
[shape,naturalDerivatives]=shapeFunctionL2(pt(1));
Xderivatives=naturalDerivatives*invJacobian;
% B matrix
B=zeros(2,2*ndof);
B(1,ndof+1:2*ndof) = Xderivatives(:)';
% K
stiffness(elementDof,elementDof)=...
stiffness(elementDof,elementDof)+...
B'*B*gaussWeights(q)*detJacobian*Ce(1,1);
force(indice+numberNodes)=force(indice+numberNodes)+...
shape*P(e,q,2)*detJacobian*gaussWeights(q);                   %%%%%%%%%%%VALEUR A VERIFIER
mass(indiceMass,indiceMass)=mass(indiceMass,indiceMass)+...
shape*shape'*gaussWeights(q)*I*rho*detJacobian;
end
end
% shear contribution for stiffness matrix
gaussLocations=[0.];                     %%%%%%%%%%%VALEUR A VERIFIER
gaussWeights=[2.];
for e=1:numberElements
indice=elementNodes(e,:);
elementDof=[ indice+numberNodes indice+2*numberNodes];
indiceMass=indice+numberNodes;
ndof=length(indice);
length_element=xx(indice(2))-xx(indice(1));
detJ0=length_element/2;invJ0=1/detJ0;
Ce = C{e};
for q=1:size(gaussWeights,1) ;
pt=gaussLocations(q,:);
[shape,naturalDerivatives]=shapeFunctionL2(pt(1));
Xderivatives=naturalDerivatives*invJacobian;
% B
B=zeros(2,2*ndof);
B(2,1:ndof) = Xderivatives(:)';
B(2,ndof+1:2*ndof) = shape;
% K
stiffness(elementDof,elementDof)=...
stiffness(elementDof,elementDof)+...
B'*B*gaussWeights(q)*detJacobian*Ce(2,2);
force(indice+2*numberNodes)=force(indice+2*numberNodes)+...
shape*P(e,q,3)*detJacobian*gaussWeights(q);                   %%%%%%%%%%%VALEUR A VERIFIER
mass(indiceMass,indiceMass)=mass(indice,indice)+shape*shape'*...
gaussWeights(q)*thickness*rho*detJacobian;
end
end
% compression contribution for stiffness matrix
gaussLocations=[0.];
gaussWeights=[2.];
for e=1:numberElements
indice=elementNodes(e,:);
elementDof=[ indice];
indiceMass=indice;
ndof=length(indice);
length_element=xx(indice(2))-xx(indice(1));
detJ0=length_element/2;invJ0=1/detJ0;
Ce = C{e};
for q=1:size(gaussWeights,1) ;
pt=gaussLocations(q,:);
[shape,naturalDerivatives]=shapeFunctionL2(pt(1));
Xderivatives=naturalDerivatives*invJacobian;
% B
B=zeros(2,ndof);
B(2,1:ndof) = Xderivatives(:)';
%B(2,ndof+1:2*ndof) = shape;
% K
stiffness(elementDof,elementDof)=...
stiffness(elementDof,elementDof)+...
B'*B*gaussWeights(q)*detJacobian*Ce(3,3);
force(indice)=force(indice)+...
shape*P(e,q,1)*detJacobian*gaussWeights(q);                   %%%%%%%%%%%VALEUR A VERIFIER
mass(indiceMass,indiceMass)=mass(indice,indice)+shape*shape'*...
gaussWeights(q)*thickness*rho*detJacobian;
end
end

function [shape,naturalDerivatives]=shapeFunctionL2(xi)
% shape function and derivatives for L2 elements
% shape : Shape functions
% naturalDerivatives: derivatives w.r.t. xi
% xi: natural coordinates (-1 ... +1)
shape=([1-xi,1+xi]/2)';
naturalDerivatives=[-1;1]/2;
%end % end function shapeFunctionL2



