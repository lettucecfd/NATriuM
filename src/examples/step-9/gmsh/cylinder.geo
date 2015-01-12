
////////////////////////////////////////////
// Constants
D    	= 1.0; 				// Diameter of cylinder
R     	= D / 2.0;			// Radius
R2    	= 5.0 * R;			// Radius of outer circle
SB    	= 0.1; 				// Dummy length
PI    	= 4.0 * Atan(1.0);	 	// Pi

oRX 	=  50.0; 			// Right x-coordinate of outer rect
oUY 	=  25.0;			// Upper y-coordinate of outer rect
oLX	= -19.0;			// Left x-coordinate of outer rect
oLY	= -25.0;			// Lower y-coordinate of outer rect

iRX 	=  5.0; 			// Right x-coordinate of inner rect
iUY 	=  4.0;				// Upper y-coordinate of inner rect
iLX	= -5.0;				// Left x-coordinate of inner rect
iLY	= -4.0;				// Lower y-coordinate of inner rect

alpha 	= Pi + Atan(-iLY/-iLX);		// Angle of lower left corner
beta  	= 2*Pi - Atan(-iLY/iRX);		// Angle of lower right corner
gamma 	= Atan(iUY/iRX);			// Angle of upper right corner
delta 	= Pi - Atan(iUY/-iLX);		// Angle of upper left corner



////////////////////////////////////////////
// Mesh parameters
NX	= 20;	// Number of cells on top/bottom boundary of inner rect
NYL	= 10;	// Number of cells on left ...
NYR	= 20;	//   ... and right inner boundary
NS	= 20;	// Number of shells in inner circle
RS	= 0.85;	// Thickness ratio between neighboring shells
NDL	= 8;	// Number of cells along left diagonals
NDR	= 8;	// Number of cells along right diagonals
NXB	= 5;	// Number of horizontal cells in the beginning section
NXE	= 20;	// Number of horizontal cells in the ending section
NYUL	= 8;	// Number of verical cells in the upper and lower section
RYUL    = 0.8;	// Height ratio in the upper and lower sections


////////////////////////////////////////////
// POINTS

// Center
Point(0)  = {  0.0,  0.0, 0.0, SB};

// Inner rectangle 
Point(1)  = {iLX, iLY, 0.0, SB};	// bottom left corner
Point(2)  = {iRX, iLY, 0.0, SB}; 	// bottom right corner
Point(3)  = {iRX, iUY, 0.0, SB}; 	// top right corner
Point(4)  = {iLX, iUY, 0.0, SB};	// top left corner

// Circle points
Point(5)  = {R*Cos(alpha), R*Sin(alpha), 0.0, SB};	//lower left
Point(6)  = {R*Cos( beta), R*Sin( beta), 0.0, SB};	//lower right
Point(7)  = {R*Cos(gamma), R*Sin(gamma), 0.0, SB};	//upper right
Point(8)  = {R*Cos(delta), R*Sin(delta), 0.0, SB};	//upper left

// Outer circle points
Point(9)  = {R2*Cos(alpha), R2*Sin(alpha), 0.0, SB};	//lower left
Point(10) = {R2*Cos( beta), R2*Sin( beta), 0.0, SB};	//lower right
Point(11) = {R2*Cos(gamma), R2*Sin(gamma), 0.0, SB};	//upper right
Point(12) = {R2*Cos(delta), R2*Sin(delta), 0.0, SB};	//upper left

// Outer rectangle; parameterized counterclockwise
Point(13)  = {oLX, oLY, 0.0, SB};	// bottom left corner
Point(14)  = {iLX+5.0, oLY, 0.0, SB};	// bottom half left
Point(15)  = {iRX+10.0, oLY, 0.0, SB};	// bottom half right
Point(16)  = {oRX, oLY, 0.0, SB}; 	// bottom right corner
Point(17)  = {oRX, iLY, 0.0, SB};	// half lower right
Point(18)  = {oRX, iUY, 0.0, SB};	// half upper right
Point(19)  = {oRX, oUY, 0.0, SB}; 	// top right corner
Point(20)  = {iRX+10.0, oUY, 0.0, SB};	// half lower right
Point(21)  = {iLX+5.0, oUY, 0.0, SB};	// half lower right
Point(22)  = {oLX, oUY, 0.0, SB};	// top left corner
Point(23)  = {oLX, iUY, 0.0, SB};	// half upper right
Point(24)  = {oLX, iLY, 0.0, SB};	// half lower right

/////////////////////////////////////////////
// LINES

// Inner rectangle
Line(1)  = { 1, 2};	//bottom
Line(2)  = { 2, 3};	//right
Line(3)  = { 3, 4};	//top
Line(4)  = { 4, 1};	//left

// Inner circle boundary
Circle(5)  = { 5, 0, 6};	//bottom
Circle(6)  = { 6, 0, 7};	//right
Circle(7)  = { 7, 0, 8};	//top
Circle(8)  = { 8, 0, 5};	//left

// Outer circle boundary
Circle(9)  = { 9, 0,10};	//bottom
Circle(10) = {10, 0,11};	//right
Circle(11) = {11, 0,12};	//top
Circle(12) = {12, 0, 9};	//left

// Diagonal lines, all directed towards center
Line(13) = { 1, 9};	//lower left, outer
Line(14) = { 9, 5};	//lower left, inner
Line(15) = { 2,10};	//lower right, outer
Line(16) = {10, 6};	//lower right, inner
Line(17) = { 3,11};	//top right, outer
Line(18) = {11, 7};	//top right, inner
Line(19) = { 4,12};	//top left, outer
Line(20) = {12, 8};	//top left, inner

// Lower left rectangle
Line(21) = {13,14};	//bottom
Line(22) = {14, 1};	//right
Line(23) = { 1,24};	//top
Line(24) = {24,13};	//left

// Lower rectangle
Line(25) = {14,15};	//bottom
Line(26) = {15, 2};	//right
//Line(27) = { 2, 1};	//top
//Line(28) = { 1,14};	//left

// Lower right rectangle
Line(29) = {15,16};	//bottom
Line(30) = {16,17};	//right
Line(31) = {17, 2};	//top
//Line(32) = { 2,15};	//left

// Right rectangle
//Line(33) = { 2,17};	//bottom
Line(34) = {17,18};	//right
Line(35) = {18, 3};	//top
//Line(36) = { 3, 2};	//left

// Upper left rectangle
//Line(37) = { 3,18};	//bottom
Line(38) = {18,19};	//right
Line(39) = {19,20};	//top
Line(40) = {20, 3};	//left

// Upper rectangle
//Line(41) = { 4, 3};	//bottom
//Line(42) = { 3,20};	//right
Line(43) = {20,21};	//top
Line(44) = {21, 4};	//left

// Upper left rectangle
Line(45) = {23, 4};	//bottom
//Line(46) = { 4,21};	//right
Line(47) = {21,22};	//top
Line(48) = {22,23};	//left

// Left rectangle
//Line(49) = {24, 1};	//bottom
//Line(50) = { 1, 4};	//right
//Line(51) = { 4,23};	//top
Line(52) = {23,24};	//left

/////////////////////////////////////////////
// Discretizations of lines

// Inner rectangle
Transfinite Line(1)  = NX;	//bottom
Transfinite Line(2)  = NYR;	//right
Transfinite Line(3)  = NX;	//top
Transfinite Line(4)  = NYL;	//left

// Inner circle boundary
Transfinite Line(5)  = NX;	//bottom
Transfinite Line(6)  = NYR;	//right
Transfinite Line(7)  = NX;	//top
Transfinite Line(8)  = NYL;	//left

// Outer circle boundary
Transfinite Line(9)  = NX;	//bottom
Transfinite Line(10) = NYR;	//right
Transfinite Line(11) = NX;	//top
Transfinite Line(12) = NYL;	//left

// Diagonal lines, all directed towards center
Transfinite Line(13) = NDL;			//lower left, outer
Transfinite Line(14) = NS Using Progression RS;	//lower left, inner
Transfinite Line(15) = NDR;           		//lower right, outer
Transfinite Line(16) = NS Using Progression RS;	//lower right, inner
Transfinite Line(17) = NDR;			//top right, outer
Transfinite Line(18) = NS Using Progression RS;	//top right, inner
Transfinite Line(19) = NDL;			//top left, outer
Transfinite Line(20) = NS Using Progression RS;	//top left, inner

// Lower left rectangle
Transfinite Line(21) = NXB;				//bottom
Transfinite Line(22) = NYUL Using Progression RYUL;	//right
Transfinite Line(23) = NXB;				//top
Transfinite Line(24) = NYUL Using Progression 1./RYUL;	//left

// Lower rectangle
Transfinite Line(25) = NX;				//bottom
Transfinite Line(26) = NYUL Using Progression RYUL;	//right
//Transfinite Line(27) = NX;				//top
//Transfinite Line(28) = NYUL;				//left

// Lower right rectangle
Transfinite Line(29) = NXE;				//bottom
Transfinite Line(30) = NYUL Using Progression RYUL;	//right
Transfinite Line(31) = NXE;				//top
//Transfinite Line(32) = NYUL;				//left

// Right rectangle
//Transfinite Line(33) = NXE;				//bottom
Transfinite Line(34) = NYR;				//right
Transfinite Line(35) = NXE;				//top
//Transfinite Line(36) = NYR;				//left

// Upper right rectangle
//Transfinite Line(37) = NXE;				//bottom
Transfinite Line(38) = NYUL Using Progression 1./RYUL;	//right
Transfinite Line(39) = NXE;				//top
Transfinite Line(40) = NYUL Using Progression RYUL;	//left

// Upper rectangle
//Transfinite Line(41) = NX;				//bottom
//Transfinite Line(42) = NYUL;				//right
Transfinite Line(43) = NX;				//top
Transfinite Line(44) = NYUL Using Progression RYUL;	//left

// Upper left rectangle
Transfinite Line(45) = NXB;				//bottom
//Transfinite Line(46) = NYUL;				//right
Transfinite Line(47) = NXB;				//top
Transfinite Line(48) = NYUL Using Progression RYUL;	//left

// Left rectangle
//Transfinite Line(49) = NXB;				//bottom
//Transfinite Line(50) = NYL;				//right
//Transfinite Line(51) = NXB;				//top
Transfinite Line(52) = NYL;				//left


/////////////////////////////////////////////
// LINE LOOPS

// Inner rectangle
Line Loop(1)  = {-13,  1, 15,  -9};	//bottom
Line Loop(2)  = {-15,  2, 17, -10};	//right
Line Loop(3)  = {-17,  3, 19, -11};	//top
Line Loop(4)  = {-19,  4, 13, -12};	//left

// Circle
Line Loop(5)  = {-14,  9, 16,  -5};	//bottom
Line Loop(6)  = {-16, 10, 18,  -6};	//right
Line Loop(7)  = {-18, 11, 20,  -7};	//top
Line Loop(8)  = {-20, 12, 14,  -8};	//left

// Outer rectangles
Line Loop(9)  = { 21, 22, 23, 24};	//lower left
Line Loop(10) = { 25, 26, -1,-22};	//lower
Line Loop(11) = { 29, 30, 31,-26};	//lower right
Line Loop(12) = {-31, 34, 35, -2};	//right
Line Loop(13) = {-35, 38, 39, 40};	//upper right
Line Loop(14) = { -3,-40, 43, 44};	//upper
Line Loop(15) = { 45,-44, 47, 48};	//upper left
Line Loop(16) = {-23, -4,-45, 52};	//left


/////////////////////////////////////////////
// Surfaces

// Inner rectangle
Ruled Surface(1)  = {1};	//bottom
Ruled Surface(2)  = {2};	//right
Ruled Surface(3)  = {3};	//top
Ruled Surface(4)  = {4};	//left

// Circle
Ruled Surface(5)  = {5};	//bottom
Ruled Surface(6)  = {6};	//right
Ruled Surface(7)  = {7};	//top
Ruled Surface(8)  = {8};	//left

// Outer rectangles
Ruled Surface(9)  = {9};	//lower left
Ruled Surface(10) = {10};	//lower
Ruled Surface(11) = {11};	//lower right
Ruled Surface(12) = {12};	//right
Ruled Surface(13) = {13};	//upper right
Ruled Surface(14) = {14};	//upper
Ruled Surface(15) = {15};	//upper left
Ruled Surface(16) = {16};	//left


/////////////////////////////////////////////
// Transfinite Surfaces -> Tria meshing

// Inner rectangle
Transfinite Surface(1)  = { 9, 1, 2,10};	//bottom
Transfinite Surface(2)  = {10, 2, 3,11};	//right
Transfinite Surface(3)  = {11, 3, 4,12};	//top
Transfinite Surface(4)  = {12, 4, 1, 9};	//left

// Circle
Transfinite Surface(5)  = { 5, 9,10, 6};	//bottom
Transfinite Surface(6)  = { 6,10,11, 7};	//right
Transfinite Surface(7)  = { 7,11,12, 8};	//top
Transfinite Surface(8)  = { 8,12, 9, 5};	//left

// Outer rectangles
Transfinite Surface(9)  = {13,14, 1,24};	//lower left
Transfinite Surface(10) = {14,15, 2, 1};	//lower
Transfinite Surface(11) = {15,16,17, 2};	//lower right
Transfinite Surface(12) = { 2,17,18, 3};	//right
Transfinite Surface(13) = { 3,18,19,20};	//upper right
Transfinite Surface(14) = { 4, 3,20,21};	//upper
Transfinite Surface(15) = {23, 4,21,22};	//upper left
Transfinite Surface(16) = {24, 1, 4,23};	//left

/////////////////////////////////////////////
// Recombine Surfaces -> Quad meshing

// Inner rectangle
Recombine Surface(1);	//bottom
Recombine Surface(2);	//right
Recombine Surface(3);	//top
Recombine Surface(4);	//left

// Inner boundary
Recombine Surface(5);	//bottom
Recombine Surface(6);	//right
Recombine Surface(7);	//top
Recombine Surface(8);	//left

// Outer rectangles
Recombine Surface(9);	//lower left
Recombine Surface(10);	//lower
Recombine Surface(11);	//lower right
Recombine Surface(12);	//right
Recombine Surface(13);	//upper right
Recombine Surface(14);	//upper
Recombine Surface(15);	//upper left
Recombine Surface(16);	//left

/////////////////////////////////////////////
// Physical IDs for reading in into deal
Physical Line(0) = {48,52,24};	//Inlet
Physical Line(1) = {30,34,38}; 	//Outlet
Physical Line(2) = {39,43,47};	//Top Wall
Physical Line(3) = {21,25,29};	//Bottom Wall
Physical Line(4) = {5,6,7,8};	//Cylinder

Physical Surface(1) = {1,2,3,4};
Physical Surface(6) = {5,6,7,8};
//Physical Surface(6) = {6};
//Phyiscal Surface(7) = {7};
//Physical Surface(8) = {8};
//Phyiscal Surface(9) = {9};
Physical Surface(10) = {9,10,11,12,13,14,15,16};
//Phyiscal Surface(11) = {11};
//Physical Surface(12) = {12};
//Phyiscal Surface(13) = {13};
//Physical Surface(14) = {14};
//Phyiscal Surface(15) = {15};
//Physical Surface(16) = {16};



