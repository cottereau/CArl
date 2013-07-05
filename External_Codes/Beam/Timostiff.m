function [x,y,K,z,F] = Timostiff(model)
%................................................................
% MATLAB codes for Finite Element Analysis
% problem16.m
% Timoshenko beam in bending
% antonio ferreira 2008
% clear memory
% E; modulus of elasticity
% G; shear modulus
% I: second moments of area
% L: length of beam
% thickness: thickness of beam
% dof = [U1 U2 ... UN V1 V2 ... VN THETA1 THETA2 ... THETAN]
% force only shear pressure
E=model.property; poisson = model.HomeFE.poisson;L = model.HomeFE.mesh.X(end,1) - model.HomeFE.mesh.X(1,1);thickness=model.HomeFE.L;
I=thickness^3/12;
%EI=E*I;
kapa=model.HomeFE.kappa;
%
P = model.load; % uniform pressure
%BC = model.HomeFE.BC;
% constitutive matrix
G=E/2/(1+poisson);
% mesh
numberElements = size(model.HomeFE.mesh.T,1);
nodeCoordinates=model.HomeFE.mesh.X;%linspace(0,L,numberElements+1);
xx=nodeCoordinates;
elementNodes=model.HomeFE.mesh.T;
% generation of coordinates and connectivities
numberNodes=size(xx,1);
% GDof: global number of degrees of freedom
GDof=3*numberNodes;
for ijk = 1:numberElements
    C{ijk} = [ E(ijk,1)*I 0 0; 0 kapa*thickness*G(ijk,1) 0;0 0 thickness*E(ijk,1)];
end
% computation of the system stiffness matrix
[stiffness,force,~]=...
    formStiffnessMassTimoshenkoBeam(GDof,numberElements,...
    elementNodes,numberNodes,xx,C,P,1,I,thickness);
% boundary conditions (simply-supported at both bords)
%fixedNodeW =[1 ; numberNodes];
%fixedNodeTX=[];
% boundary conditions (clamped at both bords)
%K = sparse(stiffness);
[x,y]=find(stiffness);
indtmp = find(stiffness);
K=stiffness(indtmp);
%F = sparse(force);
z = find(force);
indtmp = find(force);
F=force(indtmp);
if isfield(model.HomeFE,'BC')&&~isempty( model.HomeFE.BC )
    
    % Dirichlet Boundary Conditions
    ind = find( model.HomeFE.BC.type == 'U' );
    Nbc = length( ind );
    Nx = max(x);
    Nf = 1;
    x = [ x ; Nx+(1:3*Nbc)'; model.HomeFE.BC.nodes(ind)';model.HomeFE.BC.nodes(ind)'+numberNodes;model.HomeFE.BC.nodes(ind)'+2*numberNodes ];
    y = [ y ; model.HomeFE.BC.nodes(ind)';model.HomeFE.BC.nodes(ind)'+numberNodes;model.HomeFE.BC.nodes(ind)'+2*numberNodes; Nx+(1:3*Nbc)' ];
    K = [ K ; ones( 6*Nbc, 1 ) ];
    z = [ z; reshape(repmat(Nx+(1:3*Nbc)',[1 Nf]),3*Nbc*Nf,1) ];
    F = [ F; reshape(repmat(model.HomeFE.BC.value(ind,:),[1 Nf]),3*Nbc*Nf,1) ];
    
    % Neumann Boundary Conditions
    ind = find( model.HomeFE.BC.type == 'F' );
    Nbc = length( ind );
    z = [ z; reshape(repmat(model.HomeFE.BC.nodes(ind)',[1 Nf]),Nbc*Nf,1) ];
    F = [ F; reshape(repmat(model.HomeFE.BC.value(ind,1),[1 Nf]),Nbc*Nf,1) ];
    
    %ind = find( model.HomeFE.BC.type == 'Fy' );
    Nbc = length( ind );
    z = [ z; reshape(repmat(model.HomeFE.BC.nodes(ind)'+numberNodes,[1 Nf]),Nbc*Nf,1) ];
    F = [ F; reshape(repmat(model.HomeFE.BC.value(ind,2),[1 Nf]),Nbc*Nf,1) ];
    
    %ind = find( model.HomeFE.BC.type == 'M' );
    Nbc = length( ind );
    z = [ z; reshape(repmat(model.HomeFE.BC.nodes(ind)'+2*numberNodes,[1 Nf]),Nbc*Nf,1) ];
    F = [ F; reshape(repmat(model.HomeFE.BC.value(ind,3),[1 Nf]),Nbc*Nf,1) ];
    
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
    indiceMass=indice+numberNodes;
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
        mass(indice,indice)=mass(indice,indice)+shape*shape'*...
            gaussWeights(q)*thickness*rho*detJacobian;
    end
end
% shear contribution for stiffness matrix
gaussLocations=[0.];                     %%%%%%%%%%%VALEUR A VERIFIER
gaussWeights=[2.];
for e=1:numberElements
    indice=elementNodes(e,:);
    elementDof=[ indice+numberNodes indice+2*numberNodes];
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
    end
end
% compression contribution for stiffness matrix
gaussLocations=[0.];
gaussWeights=[2.];
for e=1:numberElements
    indice=elementNodes(e,:);
    elementDof=[ indice];
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



