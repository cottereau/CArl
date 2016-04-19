// Gmsh project created on Wed Mar 23 16:39:34 2016

w = 0.1;

Point(1) = { 0, 0, 0, w};
Point(2) = { 0, 2, 0, w};
Point(3) = { 2, 2, 0, w};
Point(4) = { 2, 1, 0, w};
Point(5) = { 1, 1, 0, w};
Point(6) = { 1, 0, 0, w};

Point(7) = { 0, 0, 1, w};
Point(8) = { 0, 2, 1, w};
Point(9) = { 2, 2, 1, w};
Point(10) = { 2, 1, 1, w};
Point(11) = { 1, 1, 1, w};
Point(12) = { 1, 0, 1, w};

Line(1) = {8, 7};
Line(2) = {7, 12};
Line(3) = {12, 11};
Line(4) = {11, 10};
Line(5) = {10, 9};
Line(6) = {9, 8};
Line(7) = {1, 6};
Line(8) = {6, 5};
Line(9) = {5, 4};
Line(10) = {4, 3};
Line(11) = {3, 2};
Line(12) = {2, 1};
Line(13) = {1, 7};
Line(14) = {2, 8};
Line(15) = {3, 9};
Line(16) = {4, 10};
Line(17) = {5, 11};
Line(18) = {6, 12};
Line Loop(19) = {1, 2, 3, 4, 5, 6};
Plane Surface(20) = {19};
Line Loop(21) = {12, 7, 8, 9, 10, 11};
Plane Surface(22) = {21};
Line Loop(23) = {14, -6, -15, 11};
Plane Surface(24) = {23};
Line Loop(25) = {5, -15, -10, 16};
Plane Surface(26) = {25};
Line Loop(27) = {9, 16, -4, -17};
Plane Surface(28) = {27};
Line Loop(29) = {8, 17, -3, -18};
Plane Surface(30) = {29};
Line Loop(31) = {7, 18, -2, -13};
Plane Surface(32) = {31};
Line Loop(33) = {12, 13, -1, -14};
Plane Surface(34) = {33};
Surface Loop(35) = {20, 34, 22, 32, 30, 28, 26, 24};
Volume(36) = {35};
Physical Volume(37) = {36};
Physical Surface(38) = {20, 34, 22, 26, 28, 30, 32, 24};
