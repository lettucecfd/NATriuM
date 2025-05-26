// Gmsh project created on Wed Oct 27 09:49:22 2021
SetFactory("OpenCASCADE");

size = 0.3;
airfoil_mesh_size=0.01;
spanwise = 8;
length = 6;

distanceFromInlet = 2;


chordLength = 1;
AoA = 3.5;
AoA = -AoA/360*2*Pi;

Alpha = 3.5;
Alpha = Alpha/360*2*Pi;
ChordHalfWidth = chordLength/2*Tan(Alpha);
SideLength = chordLength/2/Cos(Alpha);

Point(1) = {-distanceFromInlet,spanwise/2,0,size};
Point(2) = {length-distanceFromInlet,spanwise/2,0,size*0.5};
Point(3) = {length-distanceFromInlet,-spanwise/2.0,0,size*0.5};
Point(4) = {-distanceFromInlet,-spanwise/2.0,0,size};

Point(66) = {5,0.0,0,size*10};


Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};
Curve Loop(9) = {5,6,7,8};



Point(11) = {0,0,0,airfoil_mesh_size};
Point(12) = {SideLength*Cos(Alpha+AoA),SideLength*Sin(Alpha+AoA),0,airfoil_mesh_size};
Point(13) = {chordLength*Cos(AoA),chordLength*Sin(AoA),0,airfoil_mesh_size};
Point(14) = {SideLength*Cos(-Alpha+AoA),SideLength*Sin(-Alpha+AoA),0,airfoil_mesh_size};
Line(15) = {11,12};
Line(16) = {12,13};
Line(17) = {13,14};
Line(18) = {14,11};
Curve Loop(19) = {15,16,17,18};

Plane Surface(20) = {9,19};
Physical Line(100) = {15,16,17,18};
Physical Line(101) = {8};
Physical Line(102) = {5};
Physical Line(103) = {6};
Physical Line(104) = {7};

Recombine Surface {20};
Physical Surface(30) = {20};




