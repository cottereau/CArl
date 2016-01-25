cl__1 = 1;
L = 2.0;
w = 0.1;

Point(1) = {-L/2, -L/2, -L/2, w};
Point(2) = {L/2, -L/2, -L/2, w};
Point(3) = {-L/2, L/2, -L/2, w};
Point(4) = {-L/2, -L/2, L/2, w};
Point(5) = {-L/2, L/2, L/2, w};
Point(6) = {L/2, -L/2, L/2, w};
Point(7) = {L/2, L/2, -L/2, w};
Point(8) = {L/2, L/2, L/2, w};

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
