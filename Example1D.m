close all
clear all
% profile on

% load models
% load Tests/Coarse1D.mat
% load Tests/Medium1D.mat
load Tests/Large1D.mat
% load Tests/StochasticCoarse1D.mat
% load Tests/StochasticMedium1D.mat

% coupling{1}.mediator.support = 2;
% model{2}.property = model{2}.property(:,1:1000);

% start computation
[ sol, out ] = CArl( model, coupling, solver );

% profile viewer

% deterministic post treatment 1D
figure; plot( model{1}.mesh.X, sol{1}, 'b' );
hold on; plot( model{2}.mesh.X, sol{2}, 'r' );

% stochastic post-treatment
return
X1 = model{1}.mesh.X;
s1 = diff( sol{1} ) / mean(diff(X1)); s1 = [s1 s1]'; s1 = s1(:);
x1 = [X1(1:(end-1)) X1(2:end)]'; x1 = x1(:);
X2 = model{2}.mesh.X;
x2 = [X2(1:(end-1)) X2(2:end)]'; x2 = x2(:);
s2 = zeros( length(x2), size(out.u.MC.u,2));
s2(1:2:end,:) = diff( out.u.MC.u, 1, 1 ) / mean(diff(X2));
s2(2:2:end,:) = s2(1:2:end,:);
ms2 = mean(s2,2);
ss2 = std(s2,[],2);
pms2 = [ ms2+ss2/sqrt(1-0.9); ms2(end:-1:1)-ss2(end:-1:1)/sqrt(1-0.9) ];
figure; fill( [x2;x2(end:-1:1)], pms2, 'y' )
hold on; plot( x1, s1, 'k-', x2, ms2, 'r--', x2, s2(:,1), 'b-' );

X1 = model{1}.mesh.X;
T1 = model{1}.mesh.T;
X2 = model{2}.mesh.X;
T2 = model{2}.mesh.T;
Ne2 = size(T2,1);
Nn2 = size(X2,1);
for i1 = 1:Ne2
    X2 = [ X2
           mean(X2(T2(i1,1:2),:),1)
           mean(X2(T2(i1,2:3),:),1)
           mean(X2(T2(i1,[3 1]),:),1) ];
    T2 = [ T2
           T2(i1,1) Nn2+1 Nn2+3
           T2(i1,2) Nn2+2 Nn2+1
           T2(i1,3) Nn2+3 Nn2+2 ];
    T2(i1,:) = Nn2 + (1:3) ;
    Nn2 = Nn2 + 3;
end
for i1 = 1:Nn2
    ind = find( abs(X2(:,1)-X2(i1,1))<1e-5 & abs(X2(:,2)-X2(i1,2))<1e-5 );
    Ni = length(ind);
    for i2 = Ni:-1:2
        T2( T2==ind(i2) ) = i1;
        T2( T2>ind(i2) ) = T2( T2>ind(i2) ) - 1;
        X2 = X2( setdiff(1:Nn2,ind(i2)), : );
        Nn2 = Nn2 - 1;
    end
    if i1 >= Nn2
        break
    end
end


