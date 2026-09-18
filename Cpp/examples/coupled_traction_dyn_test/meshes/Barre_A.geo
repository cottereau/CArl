// Definition de la barre comme dans le papier de Michael Brun

L = 2.6;
l = 0.1;

lc = 0.1;

// Definition du carre de base
Point (1) = {0,0,0,lc};
Point (2) = {l,0,0,lc};
Point (3) = {l,l,0,lc};
Point (4) = {0,l,0,lc};

Line (1) = {1,2};
Line (2) = {2,3};
Line (3) = {3,4};
Line (4) = {4,1};

// Definition de la surface de base
Line Loop (1) = {1,2,3,4};
Plane Surface (1) = {1};

// Pour avoir de l'hexa
Transfinite Line{1,2,3,4} = 2;
Transfinite Surface{1};
Recombine Surface{1};

//Now make 3D by extrusion.
newEntities[] =
Extrude { 0,0,L }
{
        Surface{1};
        Layers{26};
        Recombine;
};

Transfinite Volume{newEntities[1]};
Recombine Volume{newEntities[1]};

// ======================================
// Physical groups surfaces
// ======================================

// Face du bas
Physical Surface(1) = {1};               // MIN_Z

// Faces laterales
Physical Surface(2) = {newEntities[2]}; // MIN_Y
Physical Surface(3) = {newEntities[3]}; // MAX_X
Physical Surface(4) = {newEntities[4]}; // MAX_Y
Physical Surface(5) = {newEntities[5]}; // MIN_X

// Face du haut
Physical Surface(6) = {newEntities[0]}; // MAX_Z

// ======================================
// Physical volume
// ======================================

Physical Volume(33) = {newEntities[1]};