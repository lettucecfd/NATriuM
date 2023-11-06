// Gmsh project created on Tue Aug  8 15:02:51 2023
// use ref-level 0 !
sl = 0.093;
y1 = 40;
y2 = 60;
y3 = 75;
xl = 150;
zl = 40;
nx = 76;
nz = 21;
ny1 = 17;
ny2 = 6;
ny3 = 3;

SetFactory("OpenCASCADE");
Point(1) = {xl*sl, y3*sl, zl*sl, 1.0};
Point(2) = {xl*sl, y3*sl, -zl*sl, 1.0};
Point(3) = {xl*sl, -y3*sl, -zl*sl, 1.0};
Point(4) = {xl*sl, -y3*sl, zl*sl, 1.0};
Point(5) = {-xl*sl, y3*sl, zl*sl, 1.0};
Point(6) = {-xl*sl, -y3*sl, zl*sl, 1.0};
Point(7) = {-xl*sl, -y3*sl, -zl*sl, 1.0};
Point(8) = {-xl*sl, y3*sl, -zl*sl, 1.0};
Point(9) = {-xl*sl, 0, -zl*sl, 1.0};
Point(10) = {xl*sl, 0, -zl*sl, 1.0};
Point(11) = {-xl*sl, 0, zl*sl, 1.0};
Point(12) = {xl*sl, 0, zl*sl, 1.0};
Point(13) = {-xl*sl, y1*sl, -zl*sl, 1.0};
Point(14) = {xl*sl, y1*sl, -zl*sl, 1.0};
Point(15) = {-xl*sl, y1*sl, zl*sl, 1.0};
Point(16) = {xl*sl, y1*sl, zl*sl, 1.0};
Point(17) = {-xl*sl, -y1*sl, -zl*sl, 1.0};
Point(18) = {xl*sl, -y1*sl, -zl*sl, 1.0};
Point(19) = {-xl*sl, -y1*sl, zl*sl, 1.0};
Point(20) = {xl*sl, -y1*sl, zl*sl, 1.0};
Point(21) = {-xl*sl, y2*sl, -zl*sl, 1.0};
Point(22) = {xl*sl, y2*sl, -zl*sl, 1.0};
Point(23) = {-xl*sl, y2*sl, zl*sl, 1.0};
Point(24) = {xl*sl, y2*sl, zl*sl, 1.0};
Point(25) = {-xl*sl, -y2*sl, -zl*sl, 1.0};
Point(26) = {xl*sl, -y2*sl, -zl*sl, 1.0};
Point(27) = {-xl*sl, -y2*sl, zl*sl, 1.0};
Point(28) = {xl*sl, -y2*sl, zl*sl, 1.0};
Line(1) = {6, 27};
Line(2) = {27, 19};
Line(3) = {19, 11};
Line(4) = {11, 15};
Line(5) = {15, 23};
Line(6) = {23, 5};
Line(7) = {5, 8};
Line(8) = {8, 21};
Line(9) = {21, 13};
Line(10) = {13, 9};
Line(11) = {9, 17};
Line(12) = {17, 25};
Line(13) = {25, 7};
Line(14) = {7, 6};
Line(15) = {6, 4};
Line(16) = {4, 28};
Line(17) = {28, 20};
Line(18) = {20, 12};
Line(19) = {12, 16};
Line(20) = {16, 24};
Line(21) = {24, 1};
Line(22) = {1, 2};
Line(23) = {2, 22};
Line(24) = {22, 14};
Line(25) = {14, 10};
Line(26) = {10, 18};
Line(27) = {18, 26};
Line(28) = {26, 3};
Line(29) = {3, 4};
Line(30) = {3, 7};
Line(31) = {8, 2};
Line(32) = {1, 5};
Curve Loop(1) = {31, -22, 32, 7};
Plane Surface(1) = {1};
Curve Loop(2) = {22, 23, 24, 25, 26, 27, 28, 29, 16, 17, 18, 19, 20, 21};
Plane Surface(2) = {2};
Curve Loop(3) = {30, 14, 15, -29};
Plane Surface(3) = {3};
Curve Loop(4) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14};
Plane Surface(4) = {4};
Curve Loop(5) = {32, -6, -5, -4, -3, -2, -1, 15, 16, 17, 18, 19, 20, 21};
Plane Surface(5) = {5};
Curve Loop(6) = {31, 23, 24, 25, 26, 27, 28, 30, -13, -12, -11, -10, -9, -8};
Plane Surface(6) = {6};
Transfinite Curve {:} = 2 Using Progression 1;
Transfinite Curve {10, 11, 4, 3, 19, 18, 26, 25} = ny1 Using Progression 1; // y-direction at 2*y1 initial shear layer thicknesses
Transfinite Curve {5, 9, 24, 20, 17, 27, 12, 2} = ny2 Using Progression 1; // y-direction at 2*y2 initial shear layer thicknesses
Transfinite Curve {6, 8, 23, 21, 16, 28, 13, 1} = ny3 Using Progression 1; // y-direction at 2*y3 initial shear layer thicknesses
Transfinite Curve {31, 32, 30, 15} = nx Using Progression 1; // x-direction
Transfinite Curve {22, 29, 14, 7} = nz Using Progression 1; // z-direction
Transfinite Surface {1} = {2, 1, 5, 8};
Transfinite Surface {2} = {1, 2, 3, 4};
Transfinite Surface {3} = {3, 4, 6, 7};
Transfinite Surface {4} = {5, 8, 7, 6};
Transfinite Surface {5} = {5, 6, 4, 1};
Transfinite Surface {6} = {8, 7, 3, 2};
Recombine Surface {1, 2, 3, 4, 5, 6};
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};
Transfinite Volume{1} = {8, 2, 1, 5, 7, 3, 4, 6};
