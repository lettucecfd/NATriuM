//FILE newUni 0deg

inlet_r      = 4;
inlet_front  = 3;
inlet_c      = .3;
outlet_c     = 6;
outlet_h     = inlet_r;
sponge_h     = inlet_r * 3;
size_foil    = 0.1;
size_in_out  = 1;
size_sponge  = 1;
point_id_front = 1;
point_id_top = 2;
point_id_back = 3;
point_id_bot = 4;
n_around     = 100;  // 100/n_coarsen; // 50/n_coarsen for plot
n_sphere     = 40;
n_inlet      = n_sphere;  // (n_coarsen+1)/2;  // +1 @1; +2 @3
channel_l    = 1.5;
channel_h    = inlet_r;
n_channel    = 98;
n_outlet_center = n_channel - n_around + 61;
progression_around = 1.05;
progression_sponge = 1.06;
progression_back = 1.02;
n_sponge_front= 27;
n_sponge_back= 26;


/// SPHERE
r = 0.1;
Point(1) = {0, 0, 0, size_foil};  // front
Point(2) = {r, r, 0, size_foil};  // top
Point(3) = {2*r, 0, 0, size_foil};  // back
Point(4) = {r, -r, 0, size_foil};  // bottom
Point(5) = {r, 0, 0, size_foil};  // center
Circle(1) = {1, 5, 2};
Circle(2) = {2, 5, 3};
Circle(3) = {3, 5, 4};
Circle(4) = {4, 5, 1};

Transfinite Curve {1:4} = n_sphere Using Progression 1;  // `Using Progression 1` not necessary!


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

Curve Loop (200) = {-203, 201, -1, 200}; // inlet loop top
Curve Loop (201) = {-200, -4, -202, 204}; // inlet loop bottom
Plane Surface(1) = {200};                                // inlet surface top
Plane Surface(2) = {201};                                // inlet surface bottom
Transfinite Surface {1} = {202, point_id_front, point_id_top, 200};           // inlet surface top
Transfinite Surface {2} = {202, 201, point_id_bot, point_id_front};          // inlet surface bottom


/// OUTLET
Point(210) = {outlet_c, outlet_h, 0, size_in_out};     // outlet top
Point(211) = {outlet_c, -outlet_h, 0, size_in_out};    // outlet bottom
Point(212) = {outlet_c, 0, 0, size_foil};              // outlet center
Line(210)  = {point_id_back, 212};                                 // outlet center line
Line(211)  = {210, 212};                               // outlet end top line
Line(212)  = {211, 212};                               // outlet end bottom line
Transfinite Curve {210} = n_outlet_center Using Progression progression_back;  // outlet center lines points
Transfinite Curve {-211, -212} = n_around Using Progression progression_around;              // outlet end lines points


/// CHANNEL
Line(220)  = {200, 210};                               // channel top line
Line(221)  = {201, 211};                               // channel bottom line
Transfinite Curve {220} = n_channel Using Progression 1;      // channel top points
Transfinite Curve {221} = n_channel Using Progression 1;      // channel bottom points


/// OUTLET / CHANNEL SURFACES
Curve Loop(222) = {-201, 220, 211, -210, -2}; // channel top
Curve Loop(223) = {202, -3, 210, -212, -221}; // channel bottom
Plane Surface(5) = {222};                               // channel top
Plane Surface(6) = {223};                               // channel bottom
Transfinite Surface {5} = {200, point_id_top, 212, 210}; // channel surface top
Transfinite Surface {6} = {point_id_bot, 201, 211, 212}; // channel surface bottom


/// SPONGE
Point(230) = {outlet_c, sponge_h, 0, size_sponge};      // sponge top back
Point(231) = {outlet_c, -sponge_h, 0, size_sponge};     // sponge bottom back
Line(230)  = {210, 230};                                // sponge outlet top line
Line(231)  = {211, 231};                                // sponge outlet bottom line
Line(232)  = {200, 230};                                // sponge top line
Line(233)  = {201, 231};                                // sponge bottom line
Curve Loop(230) = {232, -230, -220};               // sponge surface top
Curve Loop(231) = {221, 231, -233};               // sponge surface bottom
Transfinite Curve {232, 233} = n_sponge_front Using Progression progression_sponge;  // sponge front/diagonal lines
Transfinite Curve {230, 231} = n_sponge_back Using Progression progression_sponge;  // sponge back lines
Plane Surface(7) = {230};                               // sponge top
Plane Surface(8) = {231};                               // sponge bottom

/// MESH SIZES
Mesh.Algorithm = 6;
Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 1; // or 3; to leave no triangles

/// BOUNDARIES
Physical Curve("Inlet", 300) = {233, 204, 203, 232};
//Physical Curve("Sponge", 301) = {234, 235};
Physical Curve("Outlet", 302) = {231, 212, 211, 230};
Physical Curve("BB_BC", 303) = {1:4};
Physical Surface(236) = {1, 2, 6, 5, 7, 8};

Mesh 2;
Save "NACA0012_0deg.msh";
