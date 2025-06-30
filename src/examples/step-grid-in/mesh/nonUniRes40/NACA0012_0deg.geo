//FILE newNonUni 0deg res40

inlet_r      = 4;
inlet_front  = 3;
inlet_c      = 0.3;
outlet_c     = 6;
outlet_h     = inlet_r;
sponge_h     = inlet_r * 3;
size_foil    = 0.1;
size_in_out  = 1;
size_sponge  = 1;
point_id_top = 27;
point_id_front = 41;
point_id_bot = 55;
n_around     = 40;
n_foil       = 2;  // number of points per foil section; 2 for just the section, 3 to split once, ...
n_inlet      = point_id_front - point_id_top + 1;  // top and bottom, each; 31 default points; additional 31-1 (30 segments) points for each n_foil above 2
channel_l    = 1.5;
channel_h    = inlet_r;
n_channel    = 140;  // number of points on top and bottom
n_outlet_center = n_channel - point_id_top + 1;
progression_around = 1.05;
progression_sponge_front = 1.05;
progression_sponge_back = 1.05;
n_sponge_front= 57;
n_sponge_back= 50;

x={1.0, 0.998459, 0.993844, 0.986185, 0.975528, 0.96194, 0.945503, 0.92632, 0.904508, 0.880203, 0.853553, 0.824724, 0.793893, 0.761249, 0.726995, 0.691342, 0.654508, 0.616723, 0.578217, 0.53923, 0.5, 0.46077, 0.421783, 0.383277, 0.345492, 0.308658, 0.273005, 0.238751, 0.206107, 0.175276, 0.146447, 0.119797, 0.095492, 0.07368, 0.054497, 0.03806, 0.024472, 0.013815, 0.006156, 0.001541, 0.0, 0.001541, 0.006156, 0.013815, 0.024472, 0.03806, 0.054497, 0.07368, 0.095492, 0.119797, 0.146447, 0.175276, 0.206107, 0.238751, 0.273005, 0.308658, 0.345492, 0.383277, 0.421783, 0.46077, 0.5, 0.53923, 0.578217, 0.616723, 0.654508, 0.691342, 0.726995, 0.761249, 0.793893, 0.824724, 0.853553, 0.880203, 0.904508, 0.92632, 0.945503, 0.96194, 0.975528, 0.986185, 0.993844, 0.998459};
y={-0.0, 0.000224, 0.000891, 0.00199, 0.003501, 0.005399, 0.007651, 0.010221, 0.013071, 0.016158, 0.019438, 0.022869, 0.026405, 0.03, 0.03361, 0.037188, 0.040686, 0.044055, 0.047242, 0.050196, 0.052862, 0.055184, 0.057108, 0.058582, 0.059557, 0.059988, 0.059841, 0.059088, 0.057712, 0.055708, 0.053083, 0.049854, 0.046049, 0.041705, 0.036867, 0.03158, 0.025893, 0.019854, 0.013503, 0.006877, 0.0, -0.006877, -0.013503, -0.019854, -0.025893, -0.03158, -0.036867, -0.041705, -0.046049, -0.049854, -0.053083, -0.055708, -0.057712, -0.059088, -0.059841, -0.059988, -0.059557, -0.058582, -0.057108, -0.055184, -0.052862, -0.050196, -0.047242, -0.044055, -0.040686, -0.037188, -0.03361, -0.03, -0.026405, -0.022869, -0.019438, -0.016158, -0.013071, -0.010221, -0.007651, -0.005399, -0.003501, -0.00199, -0.000891, -0.000224};
nx = #x[];

For i In {0:nx-1}
  Point(i+1) = {x[i], y[i], 0, 1};
EndFor
For i In {1:nx-1}
  Line(i) = {i, i+1};
EndFor
Line(nx) = {nx, 1};
Transfinite Curve {1:nx} = n_foil Using Progression 1;  // `Using Progression 1` not necessary!


/// INLET
Point(200) = {0, inlet_r, 0, size_in_out};       // inlet top
Point(201) = {0, -inlet_r, 0, size_in_out};      // inlet bottom
Point(202) = {inlet_c-inlet_front, 0, 0, size_in_out};     // inlet front
Point(203) = {inlet_c, 0, 0};                          // inlet center
Line(200)  = {point_id_front, 202};                               // inlet front line
Line(201)  = {200, point_id_top};                      // inlet top line
Line(202)  = {201, point_id_bot};                      // inlet bottom line
Ellipse(203) = {200, 203, 203, 202};                   // inlet circle line top
Ellipse(204) = {202, 203, 203, 201};                   // inlet circle line bottom
Transfinite Curve {200, -201, -202} = n_around Using Progression progression_around;  // inlet lines points
Transfinite Curve {203, 204} = n_inlet Using Progression 1;        // inlet circle points
Curve Loop (200) = {-203, 201, point_id_top:point_id_front-1, 200}; // inlet loop top
Curve Loop (201) = {-200, point_id_front:point_id_bot-1, -202, -204}; // inlet loop bottom
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
Transfinite Curve {210, 213, 214} = n_outlet_center Using Progression 1;  // outlet center lines points
Transfinite Curve {-211, -212} = n_around Using Progression progression_around;              // outlet end lines points

/// CHANNEL
Line(220)  = {200, 210};                               // channel top line
Line(221)  = {201, 211};                               // channel bottom line
Transfinite Curve {220} = n_channel Using Progression 1;      // channel top points
Transfinite Curve {221} = n_channel Using Progression 1;      // channel bottom points

/// OUTLET / CHANNEL SURFACES
Curve Loop(222) = {-201, 220, 211, -210, 1:point_id_top-1}; // channel top
Curve Loop(223) = {202, point_id_bot:nx, 210, -212, -221}; // channel bottom
Plane Surface(5) = {222};                               // channel top
Plane Surface(6) = {223};                               // channel bottom
Transfinite Surface {5} = {point_id_top, 212, 210, 200}; // channel surface top
Transfinite Surface {6} = {point_id_bot, 201, 211, 212}; // channel surface bottom

/// SPONGE
Point(230) = {outlet_c, sponge_h, 0, size_sponge};      // sponge top back
Point(231) = {outlet_c, -sponge_h, 0, size_sponge};     // sponge bottom back
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
//Mesh.Algorithm = 6;
Mesh.RecombineAll = 1;
Mesh.SubdivisionAlgorithm = -1;
Mesh.RecombinationAlgorithm = 1; // or 3; to leave no triangles

/// BOUNDARIES
Physical Curve(300) = {233, 204, 203, 232};  // "Inlet", 
Physical Curve(302) = {231, 212, 211, 230};  // "Outlet", 
Physical Curve(303) = {1:nx};  // "BB_BC", 
Physical Surface(236) = {1, 2, 6, 5, 7, 8};

Mesh 2;
Save "NACA0012_0deg.msh";
