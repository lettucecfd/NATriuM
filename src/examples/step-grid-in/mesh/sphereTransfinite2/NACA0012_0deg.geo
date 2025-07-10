//FILE newUni 0deg

inlet_r      = 4;
inlet_front  = 3;
inlet_c      = .3;
outlet_c     = 6;
outlet_h     = inlet_r;
sponge_h     = inlet_r * 3;
sf    = 0.1;
size_in_out  = 1;
size_sponge  = 1;
point_id_front = 1;
point_id_top = 2;
point_id_back = 3;
point_id_bot = 4;
n_around     = 100;
n_sphere_q   = 20;
n_sphere     = 2*n_sphere_q;
n_inlet      = n_sphere;
channel_l    = 1.5;
channel_h    = inlet_r;
progression_around = 1.05;
progression_sponge = 1.06;
progression_back = 1.02;
n_sponge_front= 27;
n_sponge_back= 26;


/// SPHERE
r = 0.1;
s2 = Sqrt(2);
Point(1) = {-r, 0, 0, sf};                   // front
Point(2) = {0, r, 0, sf};                     // top
Point(3) = {+r, 0, 0, sf};                   // back
Point(4) = {0, -r, 0, sf};                    // bottom
Point(5) = {0, 0, 0, sf};                     // center
Point(6) = {r/Sqrt(2), r/Sqrt(2), 0, sf};   // back top
Point(7) = {r/Sqrt(2), -r/Sqrt(2), 0, sf};  // back bottom
Circle(1) = {1, 5, 2};
Circle(2) = {2, 5, 6};
Circle(3) = {3, 5, 7};
Circle(4) = {4, 5, 1};
Circle(5) = {6, 5, 3};
Circle(6) = {7, 5, 4};

Point(11) = {0, 1, 0, sf};   // top
Point(12) = {1, 1, 0, sf};   // back top
Point(13) = {1, -1, 0, sf};  // back bottom
Point(14) = {0, -1, 0, sf};  // bottom

Line(11) = {11,12};
Line(12) = {12,13};
Line(13) = {13,14};

Transfinite Curve {1,4} = n_sphere;
Transfinite Curve {2,3,5,6} = n_sphere_q;
Transfinite Curve {2,3,5,6} = n_sphere_q;


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
Point(212) = {outlet_c, 0, 0, sf};              // outlet center
Line(210)  = {point_id_back, 212};                                 // outlet center line
Line(211)  = {210, 212};                               // outlet end top line
Line(212)  = {211, 212};                               // outlet end bottom line
Line(213)  = {6, 210};                               // outlet diagonal top line
Line(214)  = {7, 211};                               // outlet diagonal bottom line
Transfinite Curve {210} = n_around Using Progression progression_back;     // outlet center lines points
Transfinite Curve {-211, -212} = n_sphere_q Using Progression progression_around;  // outlet end lines points
Transfinite Curve {213, 214} = n_around Using Progression progression_around;  // outlet end lines points


/// CHANNEL
Line(220)  = {200, 210};                    // channel top line
Line(221)  = {201, 211};                    // channel bottom line
Transfinite Curve {220, 221} = n_sphere_q;  // channel outlet points


/// OUTLET / CHANNEL SURFACES
Curve Loop(222) = {-201, 220, -213, -2};      // channel top
Curve Loop(223) = {202, -6, 214, -221};       // channel bottom
Curve Loop(224) = {-213, 5, 210, -211};       // channel top back
Curve Loop(225) = {214, 212, -210, 3};        // channel bottom back
Plane Surface(3) = {222};                     // channel top
Plane Surface(4) = {223};                     // channel bottom
Plane Surface(5) = {224};                     // channel top back
Plane Surface(6) = {225};                     // channel bottom back
Transfinite Surface {3} = {200, 2, 6, 210};   // channel surface top
Transfinite Surface {4} = {4, 7, 211, 201};   // channel surface bottom
Transfinite Surface {5} = {3, 6, 210, 212};   // channel surface top back
Transfinite Surface {6} = {7, 3, 212, 211};   // channel surface bottom back


/// SPONGE
Point(230) = {outlet_c, sponge_h, 0, size_sponge};      // sponge top back
Point(231) = {outlet_c, -sponge_h, 0, size_sponge};     // sponge bottom back
Line(230)  = {210, 230};                                // sponge outlet top line
Line(231)  = {211, 231};                                // sponge outlet bottom line
Line(232)  = {200, 230};                                // sponge top line
Line(233)  = {201, 231};                                // sponge bottom line
Curve Loop(230) = {232, -230, -220};                    // sponge surface top
Curve Loop(231) = {221, 231, -233};                     // sponge surface bottom
Transfinite Curve {232, 233} = n_sponge_front Using Progression progression_sponge;  // sponge front/diagonal lines
Transfinite Curve {230, 231} = n_sponge_back Using Progression progression_sponge;  // sponge back lines
Plane Surface(7) = {230};                               // sponge top
Plane Surface(8) = {231};                               // sponge bottom

/// MESH SIZES
Mesh.Algorithm = 6;
//Mesh.RecombineAll = 1;
//Mesh.RecombinationAlgorithm = 1; // or 3; to leave no triangles

/// BOUNDARIES
Physical Curve("Inlet", 300) = {233, 204, 203, 232};
Physical Curve("Outlet", 302) = {231, 212, 211, 230};
Physical Curve("BB_BC", 303) = {1:4};
Physical Surface(236) = {1, 2, 6, 5, 7, 8};

//Mesh 2;
//Save "NACA0012_0deg.msh";
