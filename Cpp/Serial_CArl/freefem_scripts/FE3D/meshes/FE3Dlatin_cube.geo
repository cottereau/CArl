// Gmsh project created on Mon Nov 23 15:03:30 2015
Point(1) = {-1.5, -1.5, -1.5, 1.0};
Point(2) = {1.5, -1.5, -1.5, 1.0};
Point(3) = {-1.5, 1.5, -1.5, 1.0};
Point(4) = {-1.5, -1.5, 1.5, 1.0};
Point(5) = {-1.5, 1.5, 1.5, 1.0};
Point(6) = {1.5, -1.5, 1.5, 1.0};
Point(7) = {1.5, 1.5, -1.5, 1.0};
Point(8) = {1.5, 1.5, 1.5, 1.0};

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

Line Loop(13) = {2, 8, 9, 7};
Plane Surface(14) = {13};
Line Loop(15) = {1, -7, -6, -5};
Plane Surface(16) = {15};
Line Loop(17) = {4, 5, -11, -10};
Plane Surface(18) = {17};
Line Loop(19) = {8, 12, -10, -3};
Plane Surface(20) = {19};
Line Loop(21) = {12, 11, 6, -9};
Plane Surface(22) = {21};
Line Loop(23) = {1, 2, 3, 4};
Plane Surface(24) = {23};

Physical Surface(25) = {14};
Physical Surface(26) = {20};
Physical Surface(27) = {18};
Physical Surface(28) = {16};
Physical Surface(29) = {22};
Physical Surface(30) = {24};

Surface Loop(31) = {14, 24, 16, 22, 20, 18};
Volume(32) = {31};
Physical Volume(33) = {32};
Transfinite Line {2, 7, 9, 8, 6, 5, 1, 3, 12, 10, 11, 4} = 10 Using Progression 1;
