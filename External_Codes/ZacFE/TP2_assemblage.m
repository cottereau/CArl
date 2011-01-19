close all;
clear all;
clc;

%---------------------------------------------------------------------
% Donn�es du probl�me
%---------------------------------------------------------------------

% Param�tres
h = 10. ;
e = 1. ;
% Materiau
Young = 210000. ;
nu = 0.3 ;
lambda = nu*Young/((1+nu)*(1-2*nu)) ;
mu = Young / (2*(1+nu)) ;
C = [lambda+2*mu lambda 0 ; lambda lambda+2*mu 0 ; 0 0 2*mu ] ;
% Effort
F = 1000. ;

% G�om�trie

% Noeuds 
Nodes{1}.V = [ 0 ; 0 ] ;
Nodes{2}.V = [ h ; 0 ] ;
Nodes{3}.V = [ 0 ; h ] ;
Nodes{4}.V = [ h ; h ] ;
Nodes{5}.V = [ 2*h ; h ] ;
Nodes{6}.V = [ 0 ; 2*h ] ;
Nodes{7}.V = [ h ; 2*h ] ;
Nodes{8}.V = [ 2*h ; 2*h ] ;

% Visualisation

figure

for k=1 : 4
    hold on
     plot([Nodes{k}.V(1,1) ; Nodes{k+1}.V(1,1)],[Nodes{k}.V(2,1) , Nodes{k+1}.V(2,1)]);
    axis equal
end
for k=6 : 7
     plot([Nodes{k}.V(1,1) ; Nodes{k+1}.V(1,1)],[Nodes{k}.V(2,1) , Nodes{k+1}.V(2,1)]);
end
for k=1 : 2
     plot([Nodes{k}.V(1,1) ; Nodes{k+2}.V(1,1)],[Nodes{k}.V(2,1) , Nodes{k+2}.V(2,1)]);
     plot([Nodes{k+3}.V(1,1) ; Nodes{k+5}.V(1,1)],[Nodes{k+3}.V(2,1) , Nodes{k+5}.V(2,1)]);
end
for k=3 : 5
     plot([Nodes{k}.V(1,1) ; Nodes{k+3}.V(1,1)],[Nodes{k}.V(2,1) , Nodes{k+3}.V(2,1)]);
end

% Elements
Elem{1}.connect = [1 2 3] ;
Elem{2}.connect = [4 3 2] ;
Elem{3}.connect = [3 4 6] ;
Elem{4}.connect = [7 6 4] ;
Elem{5}.connect = [4 5 7] ;
Elem{6}.connect = [8 7 5] ;

% Matrices B ;
Bel = [-1/h 0 1/h 0 0 0; 0 -1/h 0 0 0 1/h; -1/(h*sqrt(2)) -1/(h*sqrt(2)) 0 1/(h*sqrt(2)) 1/(h*sqrt(2)) 0];
for i=1:3
    B{2*i-1} = Bel ;
    B{2*i} = - Bel ;
end

% Matrice de rigidite elementaire
Ke = (e/2)*h^2*Bel'*C*Bel ;
keyboard

% Assemblage
Nddl = 2*length(Nodes);
K = zeros(Nddl) ;
for j = 1:length(Elem)
    el1 = Elem{j}.connect(1) ;
    el2 = Elem{j}.connect(2) ;
    el3 = Elem{j}.connect(3) ;
    K(2*el1-1:2*el1,2*el1-1:2*el1) = K(2*el1-1:2*el1,2*el1-1:2*el1)+ Ke(1:2,1:2) ;
    K(2*el2-1:2*el2,2*el2-1:2*el2) = K(2*el2-1:2*el2,2*el2-1:2*el2)+ Ke(3:4,3:4) ;
    K(2*el3-1:2*el3,2*el3-1:2*el3) = K(2*el3-1:2*el3,2*el3-1:2*el3)+ Ke(5:6,5:6) ;
    K(2*el1-1:2*el1,2*el2-1:2*el2) = K(2*el1-1:2*el1,2*el2-1:2*el2)+ Ke(1:2,3:4) ;
    K(2*el1-1:2*el1,2*el3-1:2*el3) = K(2*el1-1:2*el1,2*el3-1:2*el3)+ Ke(1:2,5:6) ;
    K(2*el2-1:2*el2,2*el3-1:2*el3) = K(2*el2-1:2*el2,2*el3-1:2*el3)+ Ke(3:4,5:6) ;
    K(2*el2-1:2*el2,2*el1-1:2*el1) = K(2*el2-1:2*el2,2*el1-1:2*el1)+ Ke(3:4,1:2) ;
    K(2*el3-1:2*el3,2*el1-1:2*el1) = K(2*el3-1:2*el3,2*el1-1:2*el1)+ Ke(5:6,1:2) ;
    K(2*el3-1:2*el3,2*el2-1:2*el2) = K(2*el3-1:2*el3,2*el2-1:2*el2)+ Ke(5:6,3:4) ;
end

% Vecteur force generalis�e
Fg = [0;0;0;0;0;0;0;0;0;e*h*F/2;0;0;0;0;0;e*h*F/2] ;

% Prise en compte des conditions d'encastrement

Kf = zeros(Nddl-4) ;
Kf(:,:) = K(5:Nddl,5:Nddl) ;
Fgf = zeros(Nddl-4,1) ;
Fgf(:) = Fg(5:Nddl) ;

% Calcul du d�placement
Uf = Kf\Fgf ;
U = zeros(Nddl,1) ;
U(5:Nddl)= Uf ;

%---------------------------------------------------------------------
% Visualisation graphique (configuration initale et d�form�e)
%---------------------------------------------------------------------

figure

for k=1 : 4
    hold on
     plot([Nodes{k}.V(1,1) ; Nodes{k+1}.V(1,1)],[Nodes{k}.V(2,1) , Nodes{k+1}.V(2,1)]);
    axis equal
end
for k=6 : 7
     plot([Nodes{k}.V(1,1) ; Nodes{k+1}.V(1,1)],[Nodes{k}.V(2,1) , Nodes{k+1}.V(2,1)]);
end
for k=1 : 2
     plot([Nodes{k}.V(1,1) ; Nodes{k+2}.V(1,1)],[Nodes{k}.V(2,1) , Nodes{k+2}.V(2,1)]);
     plot([Nodes{k+3}.V(1,1) ; Nodes{k+5}.V(1,1)],[Nodes{k+3}.V(2,1) , Nodes{k+5}.V(2,1)]);
end
for k=3 : 5
     plot([Nodes{k}.V(1,1) ; Nodes{k+3}.V(1,1)],[Nodes{k}.V(2,1) , Nodes{k+3}.V(2,1)]);
end

for k=1:length(Nodes)
    Nodes{k}.V(1,1) =  Nodes{k}.V(1,1) + U(2*k-1);
    Nodes{k}.V(2,1) =  Nodes{k}.V(2,1) + U(2*k);
end

for k=1 : 4
     plot([Nodes{k}.V(1,1) ; Nodes{k+1}.V(1,1)],[Nodes{k}.V(2,1) , Nodes{k+1}.V(2,1)],'r');
end
for k=6 : 7
     plot([Nodes{k}.V(1,1) ; Nodes{k+1}.V(1,1)],[Nodes{k}.V(2,1) , Nodes{k+1}.V(2,1)],'r');
end
for k=1 : 2
     plot([Nodes{k}.V(1,1) ; Nodes{k+2}.V(1,1)],[Nodes{k}.V(2,1) , Nodes{k+2}.V(2,1)],'r');
     plot([Nodes{k+3}.V(1,1) ; Nodes{k+5}.V(1,1)],[Nodes{k+3}.V(2,1) , Nodes{k+5}.V(2,1)],'r');
end
for k=3 : 5
     plot([Nodes{k}.V(1,1) ; Nodes{k+3}.V(1,1)],[Nodes{k}.V(2,1) , Nodes{k+3}.V(2,1)],'r');
end

%---------------------------------------------------------------------
% Post-traitement (sigma Von Mises)
%---------------------------------------------------------------------

Svmt = zeros(length(Elem),1) ;
figure
hold on
axis equal
colorbar

for j = 1:length(Elem)
    el1 = Elem{j}.connect(1) ;
    el2 = Elem{j}.connect(2) ;
    el3 = Elem{j}.connect(3) ;
    Uel{j}(1:2,1) = U(2*el1-1:2*el1,1) ;
    Uel{j}(3:4,1) = U(2*el2-1:2*el2,1) ;
    Uel{j}(5:6,1) = U(2*el3-1:2*el3,1) ;
    sigm{j} = C*B{j}*Uel{j} ;
    Svmt(j) = sqrt((sigm{j}(1)-sigm{j}(2))^2 + (3/2)*sqrt(1/2)*sigm{j}(3)^2 - sqrt(1/2)*sigm{j}(3)*(sigm{j}(1)+sigm{j}(2)));
    fill([Nodes{el1}.V(1,1) ; Nodes{el2}.V(1,1) ;  Nodes{el3}.V(1,1)],[Nodes{el1}.V(2,1) ; Nodes{el2}.V(2,1) ;  Nodes{el3}.V(2,1)],Svmt(j))
end



