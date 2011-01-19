clear all;
close all;
% clc;

% node = [ 0   0
%          0   1
%          1   1             
%          1   0] ;
% hdata.hmax = 0.5 ;
% [ X1, T1 ] = mesh2d( node , [], [] ) ;
% 
% n=find(X1(:,1)~=0) ;
% m=find(X1(:,1)==1) ;

load MeshArlequin_1.5_0.6_coarse05.mat
n=find(X1(:,1)~=-3) ;
m=find(X1(:,1)==3) ;

[Nn d]= size(X1) ;
[Ne nnode] = size(T1) ;

Fi = spalloc(Nn,1,0) ;
Fi(m,1) = 1 ;

sol = zeros(Nn,1) ;

% E1 = ones(Nn,1) ;

% DEFINE PARAMETER FIELDS
e1 = 1 ;
L1 = [ 2 2 ] ;
d1 = 0. ;

% STOCHASTIC FINITE ELEMENT OPTIONS
opts.CorrelationTrace = 0.9999 ;
opts.MonteCarloTrials = 1 ;

E1 = e1 * makePropertyMap( X1, L1, d1, opts );

[Ki,B,b]=genmat2(X1,T1,E1,zeros(Nn,1),zeros(Nn,1));

% keyboard

K = Ki(n,n) ;
F = Fi(n,1) ;
soli = K \ F ;

sol(n,1) = soli ;

figure
trisurf( T1, X1(:,1), X1(:,2), sol) ;
figure
trisurf( T1, X1(:,1), X1(:,2), E1) ;