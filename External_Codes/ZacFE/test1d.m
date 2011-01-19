% 
% Test 1D of the multiscale method
% Cantilever beam
%

close all
clear all
clc

tic ;

% Geometry
L = 10 ;

% Elements
Ng = 2 ;
Nfk = 5 ;
Nf = Nfk*Ng ;
h = L/Nf ;
X1 = [[0:h:L]', zeros(Nf+1,1)] ;

% for Gauss quadrature
elem = 2 ; % 1D element
ngaus = 5 ;

% DEFINE PARAMETER FIELDS
e1 = 1 ;
L1 = [ 1 1 ] ;
d1 = 0.00000002 ;

% STOCHASTIC FINITE ELEMENT OPTIONS
opts.CorrelationTrace = 0.9999 ;
opts.MonteCarloTrials = 1 ;
Nmontecarlo = opts.MonteCarloTrials ;

% Boundary condition
% force field at the right end
Frend = 1 ;
Fi = zeros(Nf+1,1) ;
Fi(Nf+1,1) = Frend ;

% left end clamped

% Boucle Monte Carlo

E1 = e1 * makePropertyMap( X1, L1, d1, opts );
Eel = [(E1(1:Nf,:)+E1(2:Nf+1,:))/2] ;

xs = [1:(Nf+1) 1:(Nf+1) 1:Nf 2:(Nf+1)] ;
ys = [1:(Nf+1) 1:(Nf+1) 2:(Nf+1) 1:Nf] ;
        
for i1 = 1:opts.MonteCarloTrials
    
    % construction of the global system
    Ks = 1/h*[ [Eel(:,i1)' 0] [0 Eel(:,i1)'] -Eel(:,i1)' -Eel(:,i1)'] ;
    Ki = sparse( xs, ys, Ks ) ;
    Bv = ones(Nf,1) ;
    Bs = h*[ [Bv'/3 0] [0 Bv'/3] Bv'/6 Bv'/6] ;
    Bi = sparse( xs, ys, Bs ) ;
    K = Bi + Ki ;
    
    [Kg,B,b]=genmat1(X1(:,1),E1',ones(1,1+Nf),zeros(1,1+Nf)) ;
    Kg = Kg+B ;
    
    % Applied BC
    K = K(2:(Nf+1),2:(Nf+1)) ;
    Kg = Kg(2:(Nf+1),2:(Nf+1)) ;
    F = Fi(2:(Nf+1),1) ;

    % Resolution
    u = K \ F ;
    ug = Kg \ F ;
    u = [0 ; u(:,1)] ;
    ug = [0 ; ug(:,1)] ;
    figure
	plot(X1(:,1),u)
	figure
	plot(X1(:,1),ug,'r')

%     % Quantity of interest
%     Uoutex(i1,1) = u(size(X1,1),1) ;

end

% figure
% plot(X1(:,1),u)
% figure
% plot(X1(:,1),E1)

T = toc 

% % Post-processing
% figure
% hist(Uoutex,100)

% save results_fin_L02_Ng3_Nfk200.mat Uoutex X1 Ng Nfk h L1 E1 Nmontecarlo T;
