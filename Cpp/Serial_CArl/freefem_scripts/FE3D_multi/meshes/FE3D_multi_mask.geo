//int BoundIDMinZ = 1;
//int BoundIDMinY = 2;
//int BoundIDMaxX = 3;
//int BoundIDMaxY = 4;
//int BoundIDMinX = 5;
//int BoundIDMaxZ = 6;

cl__1 = 1;
L1 = 4.0;
L2 = 2.0;
L3 = 1.6;
w = 0.4;

// Outer cube
Point(1) = {-L1/2, -L1/2, -L1/2, w};
Point(2) = {L1/2, -L1/2, -L1/2, w};
Point(3) = {-L1/2, L1/2, -L1/2, w};
Point(4) = {-L1/2, -L1/2, L1/2, w};
Point(5) = {-L1/2, L1/2, L1/2, w};
Point(6) = {L1/2, -L1/2, L1/2, w};
Point(7) = {L1/2, L1/2, -L1/2, w};
Point(8) = {L1/2, L1/2, L1/2, w};

// Inner cube
Point(9) = {-L2/2, -L2/2, -L2/2, w};
Point(10) = {L2/2, -L2/2, -L2/2, w};
Point(11) = {-L2/2, L2/2, -L2/2, w};
Point(12) = {-L2/2, -L2/2, L2/2, w};
Point(13) = {-L2/2, L2/2, L2/2, w};
Point(14) = {L2/2, -L2/2, L2/2, w};
Point(15) = {L2/2, L2/2, -L2/2, w};
Point(16) = {L2/2, L2/2, L2/2, w};

// Mediator cube
Point(17) = {-L3/2, -L3/2, -L3/2, w};
Point(18) = {L3/2, -L3/2, -L3/2, w};
Point(19) = {-L3/2, L3/2, -L3/2, w};
Point(20) = {-L3/2, -L3/2, L3/2, w};
Point(21) = {-L3/2, L3/2, L3/2, w};
Point(22) = {L3/2, -L3/2, L3/2, w};
Point(23) = {L3/2, L3/2, -L3/2, w};
Point(24) = {L3/2, L3/2, L3/2, w};

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
Line(25) = {13, 12};
Line(26) = {12, 14};
Line(27) = {14, 16};
Line(28) = {16, 13};
Line(29) = {16, 15};
Line(30) = {14, 10};
Line(31) = {10, 15};
Line(32) = {10, 9};
Line(33) = {9, 11};
Line(34) = {11, 15};
Line(35) = {11, 13};
Line(36) = {9, 12};
Line(37) = {20, 22};
Line(38) = {22, 18};
Line(39) = {18, 17};
Line(40) = {17, 20};
Line(41) = {20, 21};
Line(42) = {21, 24};
Line(43) = {24, 22};
Line(44) = {24, 23};
Line(45) = {23, 18};
Line(46) = {19, 23};
Line(47) = {19, 17};
Line(48) = {19, 21};
Line Loop(49) = {28, 25, 26, 27};
Plane Surface(50) = {49};
Line Loop(51) = {27, 29, -31, -30};
Plane Surface(52) = {51};
Line Loop(53) = {34, -31, 32, 33};
Plane Surface(54) = {53};
Line Loop(55) = {35, 25, -36, 33};
Plane Surface(56) = {55};
Line Loop(57) = {35, -28, 29, -34};
Plane Surface(58) = {57};
Line Loop(59) = {36, 26, 30, 32};
Plane Surface(60) = {59};
Line Loop(61) = {41, 42, 43, -37};
Plane Surface(62) = {61};
Line Loop(63) = {44, 45, -38, -43};
Plane Surface(64) = {63};
Line Loop(65) = {46, 45, 39, -47};
Plane Surface(66) = {65};
Line Loop(67) = {48, -41, -40, -47};
Plane Surface(68) = {67};
Line Loop(69) = {48, 42, 44, -46};
Plane Surface(70) = {69};
Line Loop(71) = {37, 38, 39, 40};
Plane Surface(72) = {71};

Surface Loop(73) = {24, 16, 14, 20, 22, 18};
Surface Loop(74) = {50, 58, 56, 60, 52, 54};
Surface Loop(75) = {62, 68, 70, 64, 66, 72};

Volume(76) = {73, 74};
Volume(77) = {74, 75};
Volume(78) = {75};
Physical Volume(79) = {76};
Physical Volume(80) = {77};
Physical Volume(81) = {78};
