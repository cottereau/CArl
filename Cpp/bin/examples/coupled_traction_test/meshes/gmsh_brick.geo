//Maillage simple de voie ferrée
//Hadrien PINAULT

//Definition des variables

dx = 10.;
dy = 2.;
dz = 2.;

x0 = +0.;
y0 = -1.;
z0 = -1.;

nex = 10;
ney = 4;
nez = 8;
Lc = nez;
//Points

//Ballast 
Point(1) = {x0,y0,z0,Lc};
Point(2) = {x0,y0+dy,z0,Lc};
Point(3) = {x0,y0+dy,z0+dz,Lc};
Point(4) = {x0,y0,z0+dz,Lc};
//Lignes
Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};
//Boucles et surfaces
Line Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};
Transfinite Line{1,3}=ney+1;
Transfinite Line{2,4}=nez+1;
// Definitions des coins
Transfinite Surface{1} = {1,2,3,4};
// Transformation des triangles en quadrangles
Recombine Surface{1};
Mesh.Smoothing = 200;
// Extrusions
Extrude {dx,0,0} {Surface{1}; Layers{10}; Recombine;}
// Physical Volumes
Physical Volume ( 1  ) = {1};
Coherence;

