cl__1 = 1e+22;
Point(1) = {0.2, 0.3, 0.2, 1e+22};
Point(2) = {1.3, 0.3, 0.2, 1e+22};
Point(3) = {1.3, 0.9, 0.2, 1e+22};
Point(4) = {0.2, 0.9, 0.2, 1e+22};
Point(5) = {0.2, 0.3, 0.7, 1e+22};
Point(6) = {1.3, 0.3, 0.7, 1e+22};
Point(10) = {1.3, 0.9, 0.7, 1e+22};
Point(14) = {0.2, 0.9, 0.7, 1e+22};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line(6) = {5, 6};
Line(7) = {6, 10};
Line(8) = {10, 14};
Line(9) = {14, 5};

Line(11) = {1, 5};
Line(12) = {2, 6};
Line(16) = {3, 10};
Line(20) = {4, 14};
Transfinite Line {1,3,6,8} = 18Using Progression 1;
Transfinite Line {2,4,7,9} = 12Using Progression 1;
Transfinite Line {11,12,16,20} = 10Using Progression 1;

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Line Loop(13) = {1, 12, -6, -11};
Ruled Surface(13) = {13};
Line Loop(17) = {2, 16, -7, -12};
Ruled Surface(17) = {17};
Line Loop(21) = {3, 20, -8, -16};
Ruled Surface(21) = {21};
Line Loop(25) = {4, 11, -9, -20};
Ruled Surface(25) = {25};
Line Loop(26) = {6, 7, 8, 9};
Plane Surface(26) = {26};
Surface Loop(1) = {1, 26, 13, 17, 21, 25};
Volume(1) = {1};
Physical Surface(1) = {1};
