// Definition de la barre comme dans le papier de Michael Brun

L_A = 2.4;
L_B = 2.4;
L_C = 0.2;
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
zoneA[] = Extrude {0,0,L_A}
{
        Surface{1};
        Layers{1};
        Recombine;
};

zoneC[] = Extrude {0,0,L_C}
{
        Surface{zoneA[0]};
        Layers{1};
        Recombine;
};

zoneB[] = Extrude {0,0,L_B}
{
        Surface{zoneC[0]};
        Layers{1};
        Recombine;
};

// Physical groups
Physical Volume(69) = {zoneA[1]};
Physical Volume(68) = {zoneC[1]};
Physical Volume(67) = {zoneB[1]};
