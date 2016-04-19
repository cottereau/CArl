// Gmsh project created on Wed Mar 23 16:39:34 2016

w = 1.2;

Point(1) = { 0.7, 0.5, 0.5, w};
Point(2) = { 1.5, 0.5, 1, w};
Point(3) = { 1.5, 0.5, 0, w};
Point(4) = { 1.5, 1.3, 0.5, w};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 2};
Line(5) = {3, 1};
Line(6) = {1, 4};

Line Loop(7) = {1, 2, 5};
Plane Surface(8) = {7};
Line Loop(9) = {2, 3, 4};
Plane Surface(10) = {9};
Line Loop(11) = {5, 6, -3};
Plane Surface(12) = {11};
Line Loop(13) = {4, -1, 6};
Plane Surface(14) = {13};
Surface Loop(15) = {14, 10, 8, 12};
Volume(16) = {15};
Physical Volume(17) = {16};

Transfinite Line {2} = 3 Using Progression 1;
Transfinite Line {1, 6, 5, 4, 3} = 2 Using Progression 1;
