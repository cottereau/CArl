cl__1 = 1;
Lx = 4.0;
Ly = 2.0;
Lz = 2.0;

xt = 3.0;
yt = 0.0;
zt = 0.0;

divx= 3;
divy = 2;
divz = 2;

w = 10;

Point(1) = {-Lx/2 + xt, -Ly/2 + yt, -Lz/2 + zt, w};
Point(2) = { Lx/2 + xt, -Ly/2 + yt, -Lz/2 + zt, w};
Point(3) = {-Lx/2 + xt,  Ly/2 + yt, -Lz/2 + zt, w};
Point(4) = {-Lx/2 + xt, -Ly/2 + yt,  Lz/2 + zt, w};
Point(5) = {-Lx/2 + xt,  Ly/2 + yt,  Lz/2 + zt, w};
Point(6) = { Lx/2 + xt, -Ly/2 + yt,  Lz/2 + zt, w};
Point(7) = { Lx/2 + xt,  Ly/2 + yt, -Lz/2 + zt, w};
Point(8) = { Lx/2 + xt,  Ly/2 + yt,  Lz/2 + zt, w};

Line(1) = {4, 6};
Line(2) = {6, 8};
Line(3) = {8, 5};
Line(4) = {5, 4};
Line(5) = {4, 1};
Line(6) = {1, 2};
Line(7) = {2, 6};
Line(8) = {8, 7};
Line(9) = {7, 2};
Line(10) = {5, 3};
Line(11) = {3, 1};
Line(12) = {7, 3};

Transfinite Line {8, 12, 6, 1, 3} = divx Using Progression 1;
Transfinite Line {2, 4, 9, 11} = divy Using Progression 1;
Transfinite Line {8, 7, 5, 10} = divz Using Progression 1;

Line Loop(14) = {2, 8, 9, 7};
Plane Surface(14) = {14};


Line Loop(16) = {1, -7, -6, -5};
Plane Surface(16) = {16};

Line Loop(18) = {4, 5, -11, -10};
Plane Surface(18) = {18};

Line Loop(20) = {8, 12, -10, -3};
Plane Surface(20) = {20};

Line Loop(22) = {12, 11, 6, -9};
Plane Surface(22) = {22};

Line Loop(24) = {1, 2, 3, 4};
Plane Surface(24) = {24};

Surface Loop(32) = {14, 24, 16, 22, 20, 18};
Volume(32) = {32};

//int BoundIDMinZ = 1;
//int BoundIDMinY = 2;
//int BoundIDMaxX = 3;
//int BoundIDMaxY = 4;
//int BoundIDMinX = 5;
//int BoundIDMaxZ = 6;

Physical Surface(1) = {22};
Physical Surface(2) = {16};
Physical Surface(3) = {14};
Physical Surface(4) = {20};
Physical Surface(5) = {18};
Physical Surface(6) = {24};

Physical Volume(33) = {32};
