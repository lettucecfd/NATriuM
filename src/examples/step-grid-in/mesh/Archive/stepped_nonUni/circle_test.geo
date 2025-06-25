SetFactory("OpenCASCADE");

// 1. Define the Cartesian mesh domain (e.g., a rectangle)
Lx = 10;
Ly = 10;
dx = 0.1; // mesh size

Point(1) = {0, 0, 0, dx};
Point(2) = {Lx, 0, 0, dx};
Point(3) = {Lx, Ly, 0, dx};
Point(4) = {0, Ly, 0, dx};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(10) = {1:4};
Plane Surface(10) = {10};
Transfinite Surface {10};

Mesh 2;

// 2. Define your complex shape using points/curves/line loops
// Example: a circle
r = 4;
Disk(20) = {5, 5, 0, r, r};

// 3. Intersect Cartesian grid surface with complex shape
//BooleanFragments{ Surface{10}; }{ Surface{20}; Delete; }
//BooleanIntersection{ Surface{10}; Delete; }{ Surface{20}; Delete; }
//BooleanDifference{ Surface{10}; Delete; }{ Surface{20}; Delete; }

Plugin(CutMesh).Shape = 20; // Reference to the cutting shape
Plugin(CutMesh).KeepInside = 1; // Keep only cells inside the shape
Plugin(CutMesh).Run = 1;