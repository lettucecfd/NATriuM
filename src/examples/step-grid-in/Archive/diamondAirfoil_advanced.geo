// Gmsh project created on Wed Oct 27 09:49:22 2021
//SetFactory("OpenCASCADE");

size = 1;
airfoil_mesh_size=1;
spanwise = 8;
length = 5;

distanceFromInlet = 2;


chordLength = 1;
AoA = 4;
Alpha = AoA;
AoA = -AoA/360*2*Pi;

Alpha = Alpha/360*2*Pi;
ChordHalfWidth = chordLength/2*Tan(Alpha);
SideLength = chordLength/2/Cos(Alpha);


///MESH INFO
number_y_cells=100;


///

Point(1) = {-distanceFromInlet,spanwise/2,0,size};


Point(2) = {length-distanceFromInlet,spanwise/2,0,size};
Point(3) = {length-distanceFromInlet,-spanwise/2.0,0,size};
Point(4) = {-distanceFromInlet,-spanwise/2.0,0,size};

Point(66) = {5,0.0,0,size*10};


//Line(5) = {1,2};
//Line(6) = {2,3};
//Line(7) = {3,4};
//Line(8) = {4,1};
//Curve Loop(9) = {5,6,7,8};



Point(11) = {0,0,0,airfoil_mesh_size};
Point(12) = {SideLength*Cos(Alpha+AoA),SideLength*Sin(Alpha+AoA),0,airfoil_mesh_size};
Point(13) = {chordLength*Cos(AoA),chordLength*Sin(AoA),0,airfoil_mesh_size};
Point(14) = {SideLength*Cos(-Alpha+AoA),SideLength*Sin(-Alpha+AoA),0,airfoil_mesh_size};
Line(15) = {11,12};
Line(16) = {12,13};
Line(17) = {13,14};
Line(18) = {14,11};
Curve Loop(19) = {15,16,17,18};

Point(121) = {0,spanwise/2,0,size};
Line (122) = {11,121};
Point(123) = {SideLength*Cos(Alpha+AoA),spanwise/2,0,size};
Line (124) = {12,123};
Point(125) = {chordLength*Cos(AoA),spanwise/2,0,size};
Line (126) = {13,125};

Point(127) = {0,-spanwise/2,0,size};
Line (128) = {11,127};
Point(129) = {SideLength*Cos(Alpha-AoA),-spanwise/2,0,size};
Line (130) = {14,129};
Point(131) = {chordLength*Cos(AoA),-spanwise/2,0,size};
Line (132) = {131,13};

Point(133) = {length-distanceFromInlet,chordLength*Sin(AoA),0,size};
Line (134) = {133,13};
Point(135) = {-distanceFromInlet,0,0,size};
Line (136) = {11,135};

Line (201) = {135,1};
Line (202) = {1,121};
Line (203) = {121,123};
Line (204) = {123,125};
Line (205) = {125,2};
Line (206) = {2,133};
Line (2060) = {133,3};
Line (207) = {3,131};
Line (208) = {131,129};
Line (209) = {129,127};
Line (210) = {127,4};
Line (211) = {4,135};

Curve Loop(301) = {201,202,-122,136};
Plane Surface(302) = {301};

Recombine Surface (302);
Transfinite Surface {302};

Curve Loop(303) = {203,-124,-15,122};
Plane Surface(304) = {303};

Curve Loop(351) = {-136,128,210,211};
Plane Surface(352) = {351};

Transfinite Line {-202,210,136} = 15 Using Progression 1.29;
Recombine Surface (352);
Transfinite Surface {352};





//Plane Surface(20) = {9,19};


//+
Transfinite Line {15} = 15 Using Progression 1.1;

Transfinite Line {203} = 15 Using Progression 1.1 ;
//+
Transfinite Line {201,122, 124, 126,-206} = number_y_cells Using Progression 1.01;

Transfinite Surface {304};
Recombine Surface {304};


Curve Loop(305) = {-18,130,209,-128};
Plane Surface(306) = {305};
Transfinite Line {-18, 209} = 15 Using Progression 1.1;

Transfinite Line {128,-211,130,-132,2060} = number_y_cells Using Progression 1.01;
Transfinite Surface {306};
Recombine Surface {306};

Curve Loop(307) = {-17,-132,208,-130};
Plane Surface(308) = {307};
Transfinite Line {17, 208} = 10 Using Progression 1;

Transfinite Surface {308};
Recombine Surface {308};

Curve Loop(309) = {-16,124,204,-126};
Plane Surface(310) = {309};
Transfinite Line {16, 204} = 10 Using Progression 1;

Transfinite Surface {310};
Recombine Surface {310};

Curve Loop(361) = {205,206,134,126};
Plane Surface(362) = {361};
Transfinite Line {205, 134} = 40 Using Progression 1;

Transfinite Surface {362};
Recombine Surface {362};

Curve Loop(371) = {-134,2060,207,132};
Plane Surface(372) = {371};
Transfinite Line {207, 134} = 40 Using Progression 1;

Transfinite Surface {372};
Recombine Surface {372};

Physical Surface(30) = {302,304,310,362,352,306,308,372};
Physical Line(100) = {15,16,17,18};
Physical Line(101) = {201,211};
Physical Line(102) = {202,203,204,205};
Physical Line(103) = {206,2060};
Physical Line(104) = {207,208,209,210};

//+
SetFactory("Built-in");
//+
SetFactory("OpenCASCADE");
