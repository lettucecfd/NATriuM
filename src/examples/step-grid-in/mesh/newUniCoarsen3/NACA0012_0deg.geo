//FILE newUniCoarsen3 0deg
//SetFactory("OpenCASCADE");

inlet_r      = 4;
inlet_front  = 3;
inlet_c      = 0.3;
outlet_c     = 6;
outlet_h     = inlet_r;
sponge_h     = inlet_r * 3;
size_foil    = 0.1;
size_in_out  = 1;
size_sponge  = 1;
n_coarsen    = 3;
n_points     = 197/n_coarsen;
point_id_top = 76/n_coarsen;
point_id_bot = (n_points+1) - (point_id_top-2);
point_id_front = 99/n_coarsen + 1;
n_around     = 100/n_coarsen; // 50/n_coarsen for plot
n_foil       = 2;
n_inlet      = point_id_front - point_id_top + (n_coarsen+1)/2;  // +1 @1; +2 @3
channel_l    = 1.5;
channel_h    = inlet_r;
n_channel    = 400/n_coarsen;  // 300/n_coarsen for plot
n_outlet_center = n_channel - point_id_top + (n_coarsen)/2;  // +1 @1; +1.5 @3
progression_around = 1.05;
progression_sponge_front = 1.1;
progression_sponge_back = 1.2;
n_sponge_front= 57/n_coarsen+7;
n_sponge_back= 50/n_coarsen;


/// FOIL
x={1.0, 0.989899, 0.979798, 0.969697, 0.959596, 0.949495, 0.939394, 0.929293, 0.919192, 0.909091, 0.89899, 0.888889, 0.878788, 0.868687, 0.858586, 0.848485, 0.838384, 0.828283, 0.818182, 0.808081, 0.79798, 0.787879, 0.777778, 0.767677, 0.757576, 0.747475, 0.737374, 0.727273, 0.717172, 0.707071, 0.69697, 0.686869, 0.676768, 0.666667, 0.656566, 0.646465, 0.636364, 0.626263, 0.616162, 0.606061, 0.59596, 0.585859, 0.575758, 0.565657, 0.555556, 0.545455, 0.535354, 0.525253, 0.515152, 0.505051, 0.494949, 0.484848, 0.474747, 0.464646, 0.454545, 0.444444, 0.434343, 0.424242, 0.414141, 0.40404, 0.393939, 0.383838, 0.373737, 0.363636, 0.353535, 0.343434, 0.333333, 0.323232, 0.313131, 0.30303, 0.292929, 0.282828, 0.272727, 0.262626, 0.252525, 0.242424, 0.232323, 0.222222, 0.212121, 0.20202, 0.191919, 0.181818, 0.171717, 0.161616, 0.151515, 0.141414, 0.131313, 0.121212, 0.111111, 0.10101, 0.0909091, 0.0808081, 0.0707071, 0.0606061, 0.0505051, 0.040404, 0.030303, 0.020202, 0.010101, 0.0, 0.010101, 0.020202, 0.030303, 0.040404, 0.0505051, 0.0606061, 0.0707071, 0.0808081, 0.0909091, 0.10101, 0.111111, 0.121212, 0.131313, 0.141414, 0.151515, 0.161616, 0.171717, 0.181818, 0.191919, 0.20202, 0.212121, 0.222222, 0.232323, 0.242424, 0.252525, 0.262626, 0.272727, 0.282828, 0.292929, 0.30303, 0.313131, 0.323232, 0.333333, 0.343434, 0.353535, 0.363636, 0.373737, 0.383838, 0.393939, 0.40404, 0.414141, 0.424242, 0.434343, 0.444444, 0.454545, 0.464646, 0.474747, 0.484848, 0.494949, 0.505051, 0.515152, 0.525253, 0.535354, 0.545455, 0.555556, 0.565657, 0.575758, 0.585859, 0.59596, 0.606061, 0.616162, 0.626263, 0.636364, 0.646465, 0.656566, 0.666667, 0.676768, 0.686869, 0.69697, 0.707071, 0.717172, 0.727273, 0.737374, 0.747475, 0.757576, 0.767677, 0.777778, 0.787879, 0.79798, 0.808081, 0.818182, 0.828283, 0.838384, 0.848485, 0.858586, 0.868687, 0.878788, 0.888889, 0.89899, 0.909091, 0.919192, 0.929293, 0.939394, 0.949495, 0.959596, 0.969697, 0.979798, 0.989899};
y={-1.66533e-17, 0.00145861, 0.00289836, 0.00431962, 0.00572277, 0.00710817, 0.00847614, 0.00982701, 0.0111611, 0.0124786, 0.0137799, 0.0150652, 0.0163347, 0.0175886, 0.0188271, 0.0200504, 0.0212587, 0.0224521, 0.0236306, 0.0247945, 0.0259438, 0.0270785, 0.0281986, 0.0293043, 0.0303955, 0.0314722, 0.0325343, 0.0335819, 0.0346147, 0.0356328, 0.0366359, 0.037624, 0.0385969, 0.0395544, 0.0404963, 0.0414223, 0.0423323, 0.0432259, 0.0441028, 0.0449627, 0.0458053, 0.0466302, 0.0474369, 0.0482251, 0.0489942, 0.0497439, 0.0504735, 0.0511826, 0.0518705, 0.0525368, 0.0531806, 0.0538013, 0.0543983, 0.0549708, 0.055518, 0.056039, 0.056533, 0.0569991, 0.0574363, 0.0578436, 0.0582199, 0.0585641, 0.058875, 0.0591513, 0.0593918, 0.0595951, 0.0597595, 0.0598837, 0.0599659, 0.0600043, 0.0599971, 0.0599423, 0.0598376, 0.0596807, 0.0594692, 0.0592004, 0.0588712, 0.0584786, 0.0580189, 0.0574884, 0.0568827, 0.0561972, 0.0554263, 0.0545642, 0.053604, 0.0525375, 0.0513557, 0.0500476, 0.0486001, 0.0469971, 0.0452189, 0.0432401, 0.0410275, 0.0385355, 0.0356993, 0.0324196, 0.0285303, 0.0237077, 0.0171188, 0.0, -0.0171188, -0.0237077, -0.0285303, -0.0324196, -0.0356993, -0.0385355, -0.0410275, -0.0432401, -0.0452189, -0.0469971, -0.0486001, -0.0500476, -0.0513557, -0.0525375, -0.053604, -0.0545642, -0.0554263, -0.0561972, -0.0568827, -0.0574884, -0.0580189, -0.0584786, -0.0588712, -0.0592004, -0.0594692, -0.0596807, -0.0598376, -0.0599423, -0.0599971, -0.0600043, -0.0599659, -0.0598837, -0.0597595, -0.0595951, -0.0593918, -0.0591513, -0.058875, -0.0585641, -0.0582199, -0.0578436, -0.0574363, -0.0569991, -0.056533, -0.056039, -0.055518, -0.0549708, -0.0543983, -0.0538013, -0.0531806, -0.0525368, -0.0518705, -0.0511826, -0.0504735, -0.0497439, -0.0489942, -0.0482251, -0.0474369, -0.0466302, -0.0458053, -0.0449627, -0.0441028, -0.0432259, -0.0423323, -0.0414223, -0.0404963, -0.0395544, -0.0385969, -0.037624, -0.0366359, -0.0356328, -0.0346147, -0.0335819, -0.0325343, -0.0314722, -0.0303955, -0.0293043, -0.0281986, -0.0270785, -0.0259438, -0.0247945, -0.0236306, -0.0224521, -0.0212587, -0.0200504, -0.0188271, -0.0175886, -0.0163347, -0.0150652, -0.0137799, -0.0124786, -0.0111611, -0.00982701, -0.00847614, -0.00710817, -0.00572277, -0.00431962, -0.00289836, -0.00145861};

nx = n_points;//#x[]/n_coarsen-1;


For i In {0: nx}
  Point(i+1) = {x[n_coarsen*i], y[n_coarsen*i], 0, 1};
EndFor

For i In {1: nx}
  Line(i) = {i, i+1};
EndFor
Line(nx+1) = {nx+1, 1};

Transfinite Curve {1:nx+1} = n_foil;


/// INLET
Point(200) = {0, inlet_r, 0, size_in_out};       // inlet top
Point(201) = {0, -inlet_r, 0, size_in_out};      // inlet bottom
Point(202) = {inlet_c-inlet_front, 0, 0, size_in_out};     // inlet front
Point(203) = {inlet_c, 0, 0};                          // inlet center
Line(200)  = {point_id_front, 202};                               // inlet front line
Line(201)  = {200, point_id_top};                      // inlet top line
Line(202)  = {201, point_id_bot};                      // inlet bottom line
Ellipse(203) = {200, 203, 200, 202};                   // inlet circle line top
Ellipse(204) = {201, 203, 201, 202};                   // inlet circle line bottom
Transfinite Curve {200, -201, -202} = n_around Using Progression progression_around;  // inlet lines points
Transfinite Curve {203, 204} = n_inlet Using Progression 1;        // inlet circle points

Curve Loop (200) = {-203, 201, point_id_top:point_id_front-1/n_coarsen, 200}; // inlet loop top
Curve Loop (201) = {-200, point_id_front:point_id_bot-1, -202, 204}; // inlet loop bottom
Plane Surface(1) = {200};                                // inlet surface top
Plane Surface(2) = {201};                                // inlet surface bottom
Transfinite Surface {1} = {202, point_id_front, point_id_top, 200};           // inlet surface top
Transfinite Surface {2} = {202, 201, point_id_bot, point_id_front};          // inlet surface bottom


/// OUTLET
Point(210) = {outlet_c, outlet_h, 0, size_in_out};     // outlet top
Point(211) = {outlet_c, -outlet_h, 0, size_in_out};    // outlet bottom
Point(212) = {outlet_c, 0, 0, size_foil};              // outlet center
Line(210)  = {1, 212};                                 // outlet center line
Line(211)  = {210, 212};                               // outlet end top line
Line(212)  = {211, 212};                               // outlet end bottom line
Transfinite Curve {210} = n_outlet_center Using Progression 1;  // outlet center lines points
Transfinite Curve {-211, -212} = n_around Using Progression progression_around;              // outlet end lines points


/// CHANNEL
Line(220)  = {200, 210};                               // channel top line
Line(221)  = {201, 211};                               // channel bottom line
Transfinite Curve {220} = n_channel Using Progression 1;      // channel top points
Transfinite Curve {221} = n_channel Using Progression 1;      // channel bottom points


/// OUTLET / CHANNEL SURFACES
Curve Loop(222) = {-201, 220, 211, -210, 1:(point_id_top-1)}; // channel top
Curve Loop(223) = {202, point_id_bot:nx+1, 210, -212, -221}; // channel bottom
Plane Surface(5) = {222};                               // channel top
Plane Surface(6) = {223};                               // channel bottom
Transfinite Surface {5} = {200, point_id_top, 212, 210}; // channel surface top
Transfinite Surface {6} = {point_id_bot, 201, 211, 212}; // channel surface bottom


/// SPONGE
Point(230) = {outlet_c, sponge_h, 0, size_sponge};      // sponge top back
Point(231) = {outlet_c, -sponge_h, 0, size_sponge};     // sponge bottom back
//Point(232) = {outlet_c/2, inlet_h + (inlet_h-sponge_h)/2, 0, size_sponge};    // sponge top front
//Point(233) = {outlet_c/2, -(inlet_h + (inlet_h-sponge_h)/2), 0, size_sponge};   // sponge bottom front
Line(230)  = {210, 230};                                // sponge outlet top line
Line(231)  = {211, 231};                                // sponge outlet bottom line
Line(232)  = {200, 230};                                // sponge top line
Line(233)  = {201, 231};                                // sponge bottom line
Curve Loop(230) = {232, -230, -220};               // sponge surface top
Curve Loop(231) = {221, 231, -233};               // sponge surface top
Transfinite Curve {232, 233} = n_sponge_front Using Progression progression_sponge_front;  // sponge front/diagonal lines
Transfinite Curve {230, 231} = n_sponge_back Using Progression progression_sponge_back;  // sponge back lines
Plane Surface(7) = {230};                               // sponge top
Plane Surface(8) = {231};                               // sponge bottom




/// MESH SIZES
//Mesh.ElementOrder = 1;
Mesh.Algorithm = 6;
//Mesh.SubdivisionAlgorithm = 1;  // 1 to subdivide as quadrangles
Mesh.RecombineAll = 1;
//Mesh.SubdivisionAlgorithm = -1;
Mesh.RecombinationAlgorithm = 1; // or 3; to leave no triangles

/// BOUNDARIES
Physical Curve("Inlet", 300) = {233, 204, 203, 232};  // "Inlet", 
// TODO: only inlet in channel, not sponge Physical Curve(300) = {204, 203};  // 
//Physical Curve("Sponge", 301) = {234, 235};
Physical Curve("Outlet", 302) = {231, 212, 211, 230};  // 
Physical Curve("BB_BC", 303) = {1:nx+1};  // 
Physical Surface(236) = {1, 2, 6, 5, 7, 8};

Mesh 2;

Save "NACA0012_0deg.msh";
