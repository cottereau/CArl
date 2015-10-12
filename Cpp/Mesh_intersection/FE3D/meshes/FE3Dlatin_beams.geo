Point(1) = {0, 0, 0};
Point(2) = {2, 0, 0};
Point(3) = {0, 1, 0};
Point(4) = {0, 0, 1};
Point(5) = {2, 1, 0};
Point(6) = {0, 1, 1};
Point(7) = {2, 0, 1};
Point(8) = {2, 1, 1};

Line(101) = {4, 7};
Line(102) = {7, 8};
Line(103) = {8, 6};
Line(104) = {6, 4};
Line(105) = {7, 2};
Line(106) = {2, 1};
Line(107) = {1, 4};
Line(108) = {8, 5};
Line(109) = {2, 5};
Line(110) = {5, 3};
Line(111) = {3, 6};
Line(112) = {3, 1};

Transfinite Line {102,104,105,107,108,109,111,112} = 3;
Transfinite Line {101,103,106,110} = 7;

Line Loop(314) = {102, 108, -109, -105};
Plane Surface(314) = {314};
Line Loop(316) = {103, -111, -110, -108};
Plane Surface(316) = {316};
Line Loop(318) = {104, -107, -112, 111};
Plane Surface(318) = {318};
Line Loop(320) = {112, -106, 109, 110};
Plane Surface(320) = {320};
Line Loop(322) = {107, 101, 105, 106};
Plane Surface(322) = {322};
Line Loop(324) = {101, 102, 103, 104};
Plane Surface(324) = {324};

Surface Loop(426) = {314, 316, 318, 320, 322, 324};
Volume(426) = {426};
Translate {1, 0, 0} {
  Volume{426};
}

