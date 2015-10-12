Icl__1 = 1;

Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {0, 1, 0, 1};
Point(4) = {0, 0, 1, 1};

Line(101) = {1, 2};
Line(102) = {1, 3};
Line(103) = {4, 1};
Line(104) = {2, 4};
Line(105) = {3, 4};
Line(106) = {2, 3};

Transfinite Line {101,102,103,104,105} = 2Using Progression 1;
Transfinite Line {106} = 3Using Progression 1;

Line Loop(214) = {101, 104, 103};
Plane Surface(314) = {214};

Line Loop(216) = {102, 105, 103};
Plane Surface(316) = {216};

Line Loop(218) = {101, 106, -102};
Plane Surface(318) = {218};

Line Loop(220) = {-104, 106, 105};
Plane Surface(320) = {220};

Surface Loop(426) = {314, 316, 318, 320};
Volume(526) = {426};
