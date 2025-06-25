inlet_r      = 4;
inlet_front  = 3;
inlet_c      = 0.3;
outlet_c     = 6;
outlet_h     = inlet_r;
sponge_h     = inlet_r * 3;
size_foil    = 0.1;
size_in_out  = 1;
size_sponge  = 1;
point_id_top = 70;
point_id_bot = 126;
n_around     = 100;
n_foil       = 2;  // number of points per foil section; 2 for just the section, 3 to split once, ...
n_inlet      = 31 + (n_foil-2)*30 + 1;  // top and bottom, each; 31 default points; additional 31-1 (30 segments) points for each n_foil above 2
channel_l    = 1.5;
channel_h    = inlet_r;
n_points_along_foil = 69 + (n_foil-2)*(69-1);  // 7069default points; additional 69-1 (segments) points for each n_foil above 2;
n_channel    = 400;  // number of points on top and bottom
n_outlet_center = n_channel + 1 - n_points_along_foil;
progression_around = 1.05;
progression_sponge_front = 1.05;
progression_sponge_back = 1.05;
n_sponge_front= 57;
n_sponge_back= 50;

//Warning : Start point 62 and end point 462 of GEO line 62 are closer than the geometrical tolerance, at position (0.306827, 0.0599951, 0)
//Warning : Start point 135 and end point 535 of GEO line 135 are closer than the geometrical tolerance, at position (0.292292, -0.0599951, 0)

x={1.0, 0.989899, 0.979798, 0.969697, 0.959596, 0.949495, 0.939394, 0.929293, 0.919192, 0.909091, 0.89899, 0.888889, 0.878788, 0.868687, 0.858586, 0.848485, 0.838384, 0.828283, 0.818182, 0.808081, 0.79798, 0.787879, 0.777778, 0.767677, 0.757576, 0.747475, 0.737374, 0.727273, 0.717172, 0.707071, 0.69697, 0.686869, 0.676768, 0.666667, 0.656566, 0.646465, 0.636364, 0.626263, 0.616162, 0.606061, 0.59596, 0.585859, 0.575758, 0.565657, 0.555556, 0.545455, 0.535354, 0.525253, 0.515152, 0.505051, 0.494949, 0.484848, 0.474747, 0.464646, 0.454545, 0.444444, 0.434343, 0.424242, 0.414141, 0.40404, 0.393939, 0.383838, 0.373737, 0.363636, 0.353535, 0.343434, 0.333333, 0.323232, 0.313131, 0.30303, 0.292929, 0.282828, 0.272727, 0.262626, 0.252525, 0.242424, 0.232323, 0.222222, 0.212121, 0.20202, 0.191919, 0.181818, 0.171717, 0.161616, 0.151515, 0.141414, 0.131313, 0.121212, 0.111111, 0.10101, 0.0909091, 0.0808081, 0.0707071, 0.0606061, 0.0505051, 0.040404, 0.030303, 0.020202, 0.010101, 0.0, 0.010101, 0.020202, 0.030303, 0.040404, 0.0505051, 0.0606061, 0.0707071, 0.0808081, 0.0909091, 0.10101, 0.111111, 0.121212, 0.131313, 0.141414, 0.151515, 0.161616, 0.171717, 0.181818, 0.191919, 0.20202, 0.212121, 0.222222, 0.232323, 0.242424, 0.252525, 0.262626, 0.272727, 0.282828, 0.292929, 0.30303, 0.313131, 0.323232, 0.333333, 0.343434, 0.353535, 0.363636, 0.373737, 0.383838, 0.393939, 0.40404, 0.414141, 0.424242, 0.434343, 0.444444, 0.454545, 0.464646, 0.474747, 0.484848, 0.494949, 0.505051, 0.515152, 0.525253, 0.535354, 0.545455, 0.555556, 0.565657, 0.575758, 0.585859, 0.59596, 0.606061, 0.616162, 0.626263, 0.636364, 0.646465, 0.656566, 0.666667, 0.676768, 0.686869, 0.69697, 0.707071, 0.717172, 0.727273, 0.737374, 0.747475, 0.757576, 0.767677, 0.777778, 0.787879, 0.79798, 0.808081, 0.818182, 0.828283, 0.838384, 0.848485, 0.858586, 0.868687, 0.878788, 0.888889, 0.89899, 0.909091, 0.919192, 0.929293, 0.939394, 0.949495, 0.959596, 0.969697, 0.979798, 0.989899};
y={-1.66533e-17, 0.00145861, 0.00289836, 0.00431962, 0.00572277, 0.00710817, 0.00847614, 0.00982701, 0.0111611, 0.0124786, 0.0137799, 0.0150652, 0.0163347, 0.0175886, 0.0188271, 0.0200504, 0.0212587, 0.0224521, 0.0236306, 0.0247945, 0.0259438, 0.0270785, 0.0281986, 0.0293043, 0.0303955, 0.0314722, 0.0325343, 0.0335819, 0.0346147, 0.0356328, 0.0366359, 0.037624, 0.0385969, 0.0395544, 0.0404963, 0.0414223, 0.0423323, 0.0432259, 0.0441028, 0.0449627, 0.0458053, 0.0466302, 0.0474369, 0.0482251, 0.0489942, 0.0497439, 0.0504735, 0.0511826, 0.0518705, 0.0525368, 0.0531806, 0.0538013, 0.0543983, 0.0549708, 0.055518, 0.056039, 0.056533, 0.0569991, 0.0574363, 0.0578436, 0.0582199, 0.0585641, 0.058875, 0.0591513, 0.0593918, 0.0595951, 0.0597595, 0.0598837, 0.0599659, 0.0600043, 0.0599971, 0.0599423, 0.0598376, 0.0596807, 0.0594692, 0.0592004, 0.0588712, 0.0584786, 0.0580189, 0.0574884, 0.0568827, 0.0561972, 0.0554263, 0.0545642, 0.053604, 0.0525375, 0.0513557, 0.0500476, 0.0486001, 0.0469971, 0.0452189, 0.0432401, 0.0410275, 0.0385355, 0.0356993, 0.0324196, 0.0285303, 0.0237077, 0.0171188, 0.0, -0.0171188, -0.0237077, -0.0285303, -0.0324196, -0.0356993, -0.0385355, -0.0410275, -0.0432401, -0.0452189, -0.0469971, -0.0486001, -0.0500476, -0.0513557, -0.0525375, -0.053604, -0.0545642, -0.0554263, -0.0561972, -0.0568827, -0.0574884, -0.0580189, -0.0584786, -0.0588712, -0.0592004, -0.0594692, -0.0596807, -0.0598376, -0.0599423, -0.0599971, -0.0600043, -0.0599659, -0.0598837, -0.0597595, -0.0595951, -0.0593918, -0.0591513, -0.058875, -0.0585641, -0.0582199, -0.0578436, -0.0574363, -0.0569991, -0.056533, -0.056039, -0.055518, -0.0549708, -0.0543983, -0.0538013, -0.0531806, -0.0525368, -0.0518705, -0.0511826, -0.0504735, -0.0497439, -0.0489942, -0.0482251, -0.0474369, -0.0466302, -0.0458053, -0.0449627, -0.0441028, -0.0432259, -0.0423323, -0.0414223, -0.0404963, -0.0395544, -0.0385969, -0.037624, -0.0366359, -0.0356328, -0.0346147, -0.0335819, -0.0325343, -0.0314722, -0.0303955, -0.0293043, -0.0281986, -0.0270785, -0.0259438, -0.0247945, -0.0236306, -0.0224521, -0.0212587, -0.0200504, -0.0188271, -0.0175886, -0.0163347, -0.0150652, -0.0137799, -0.0124786, -0.0111611, -0.00982701, -0.00847614, -0.00710817, -0.00572277, -0.00431962, -0.00289836, -0.00145861};

//iStepped = 400;
nx = #x[];

/// FOIL
For i In {0: nx-2}
  Point(2*i) = {x[i], y[i], 0, 1};
  Point(2*i+1) = {x[i], y[i+1], 0, 1};
EndFor

Point(2*(nx-1)) = {x[nx-1], y[nx-1], 0, 1};
Point(2*(nx-1)+1) = {x[nx-1], y[0], 0, 1};

For i In {0: nx-2}
  Line(2*i) = {2*i, 2*i+1};
  Line(2*i+1) = {2*i+1, 2*(i+1)};
EndFor

Line(2*(nx-1)) = {2*(nx-1), 2*(nx-1)+1};
Line(2*(nx-1)+1) = {2*(nx-1)+1, 0};

Transfinite Curve {0:(2*nx-1)} = n_foil Using Progression 1; // foil points


/// INLET
Point(400) = {0, inlet_r, 0, size_in_out};       // inlet top
Point(401) = {0, -inlet_r, 0, size_in_out};      // inlet bottom
Point(402) = {inlet_c-inlet_front, 0, 0, size_in_out};     // inlet front
Point(403) = {inlet_c, 0, 0};                          // inlet center
Line(400)  = {198, 402};                               // inlet front line
Line(401)  = {400, 2*point_id_top};                      // inlet top line
Line(402)  = {401, 2*point_id_bot+1};                      // inlet bottom line
Ellipse(403) = {400, 403, 403, 402};                   // inlet circle line top
Ellipse(404) = {402, 403, 403, 401};                   // inlet circle line bottom
Transfinite Curve {400, -401, -402} = n_around Using Progression progression_around;  // inlet lines points
Transfinite Curve {403} = n_inlet+1 Using Progression 1;        // inlet circle points
Transfinite Curve {404} = n_inlet Using Progression 1;        // inlet circle points
Curve Loop (500) = {-403, 401, 140:197, 400}; // inlet loop top
Curve Loop (501) = {-400, 198:252, -402, -404}; // inlet loop bottom
Plane Surface(1) = {500};                                // inlet surface top
Plane Surface(2) = {501};                                // inlet surface bottom
//Transfinite Surface {1} = {402, 100, 70, 500};           // inlet surface top
//Transfinite Surface {2} = {402, 401, 130, 100};          // inlet surface bottom

/// OUTLET
Point(410) = {outlet_c, outlet_h, 0, size_in_out};     // outlet top
Point(411) = {outlet_c, -outlet_h, 0, size_in_out};    // outlet bottom
Point(412) = {outlet_c, 0, 0, size_foil};              // outlet center
Line(410)  = {0, 412};                                 // outlet center line
Line(411)  = {410, 412};                               // outlet end top line
Line(412)  = {411, 412};                               // outlet end bottom line
Transfinite Curve {410, 413, 414} = n_outlet_center Using Progression 1;  // outlet center lines points
Transfinite Curve {-411} = n_around Using Progression progression_around;              // outlet end lines points
Transfinite Curve {-412} = n_around-1 Using Progression progression_around;              // outlet end lines points

/// CHANNEL
Line(420)  = {400, 410};                               // channel top line
Line(421)  = {401, 411};                               // channel bottom line
Transfinite Curve {420} = n_channel/8 Using Progression 1;      // channel top points
Transfinite Curve {421} = n_channel/8 Using Progression 1;      // channel bottom points


/// OUTLET / CHANNEL SURFACES
Curve Loop(522) = {-401, 420, 411, -410, 0:139}; // channel top
Curve Loop(523) = {402, (2*point_id_bot+1):(2*nx-1), 410, -412, -421}; // channel bottom
Plane Surface(5) = {522};                               // channel top
Plane Surface(6) = {523};                               // channel bottom
//Transfinite Surface {5} = {2*point_id_top, 412, 410, 400}; // channel surface top
//Transfinite Surface {6} = {2*point_id_bot+1, 401, 411, 412}; // channel surface bottom

/// SPONGE
Point(430) = {outlet_c, sponge_h, 0, size_sponge};      // sponge top back
Point(431) = {outlet_c, -sponge_h, 0, size_sponge};     // sponge bottom back
//Point(432) = {outlet_c/2, inlet_h + (inlet_h-sponge_h)/2, 0, size_sponge};    // sponge top front
//Point(433) = {outlet_c/2, -(inlet_h + (inlet_h-sponge_h)/2), 0, size_sponge};   // sponge bottom front
Line(430)  = {410, 430};                                // sponge outlet top line
Line(431)  = {411, 431};                                // sponge outlet bottom line
Line(432)  = {400, 430};                                // sponge top line
Line(433)  = {401, 431};                                // sponge bottom line
Curve Loop(530) = {432, -430, -420};               // sponge surface top
Curve Loop(531) = {421, 431, -433};               // sponge surface top
Transfinite Curve {432, 433} = n_sponge_front Using Progression progression_sponge_front;  // sponge front/diagonal lines
Transfinite Curve {430, 431} = n_sponge_back Using Progression progression_sponge_back;  // sponge back lines
Plane Surface(7) = {530};                               // sponge top
Plane Surface(8) = {531};                               // sponge bottom

//Transfinite Surface{7};
//Recombine Surface{1,2,5,6,7};

/// MESH SIZES
Mesh.ElementOrder = 1;
Mesh.Algorithm = 6;
Mesh.RecombineAll = 1;
//Mesh.SubdivisionAlgorithm = 1;
Mesh.RecombinationAlgorithm = 1; // or 3; to leave no triangles

/// BOUNDARIES
Physical Curve(300) = {433, 404, 403, 432};  // "Inlet", 
// TODO: only inlet in channel, not sponge Physical Curve(300) = {204, 203};  // "Inlet", 
//Physical Curve("Sponge", 301) = {234, 235};
Physical Curve(302) = {431, 412, 411, 430};  // "Outlet", 
Physical Curve(303) = {0,(2*nx-2)};  // "BB_BC", 
Physical Surface(236) = {1, 2, 6, 5, 7, 8};

Mesh 2;
RecombineMesh;
//RecombineMesh;
//RecombineMesh;

Save "NACA0012_0deg.msh";
