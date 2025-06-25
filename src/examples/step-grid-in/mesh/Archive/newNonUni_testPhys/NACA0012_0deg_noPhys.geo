SetFactory("OpenCASCADE");

Point(0) = {0, 0, 0, 1};
Point(1) = {0, 10, 0, 1};
Point(2) = {10, 10, 0, 1};
Point(3) = {10, 0, 0, 1};

For i In {0: 2}
  Line(i) = {i, i+1};
EndFor
Line(3) = {3,0};

Transfinite Curve {0:3} = 4 Using Progression 1;
Curve Loop(1) = {0:3};

Plane Surface(1) = {1};
Transfinite Surface{1};

/// MESH SIZES
//Mesh.ElementOrder = 1;
//Mesh.Algorithm = 6;
//Mesh.SubdivisionAlgorithm = 1;  // 1 to subdivide as quadrangles
//Mesh.RecombineAll = 1;
//Mesh.SubdivisionAlgorithm = 0;
//Mesh.RecombinationAlgorithm = 1; // or 3; to leave no triangles

/// BOUNDARIES


Mesh 2;

Save "NACA0012_0deg_noPhys.msh";

Physical Curve("Inlet", 300) = {0};  // 

Save "NACA0012_0deg.msh";