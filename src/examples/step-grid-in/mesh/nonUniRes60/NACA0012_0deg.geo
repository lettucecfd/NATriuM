//FILE newNonUni 0deg res60

inlet_r      = 4;
inlet_front  = 3;
inlet_c      = 0.3;
outlet_c     = 6;
outlet_h     = inlet_r;
sponge_h     = inlet_r * 3;
size_foil    = 0.1;
size_in_out  = 1;
size_sponge  = 1;
point_id_top = 40;
point_id_front = 60;
point_id_bot = 80;
n_around     = 60;
n_foil       = 2;  // number of points per foil section; 2 for just the section, 3 to split once, ...
n_inlet      = point_id_front - point_id_top + 1;  // top and bottom, each; 31 default points; additional 31-1 (30 segments) points for each n_foil above 2
channel_l    = 1.5;
channel_h    = inlet_r;
n_channel    = 210;  // number of points on top and bottom
n_outlet_center = n_channel - point_id_top + 1;
progression_around = 1.05;
progression_sponge_front = 1.05;
progression_sponge_back = 1.05;
n_sponge_front= 57;
n_sponge_back= 50;

x={1.0, 0.997261, 0.993844, 0.989074, 0.982963, 0.975528, 0.96679, 0.956773, 0.945503, 0.933013, 0.919335, 0.904508, 0.888573, 0.871572, 0.853553, 0.834565, 0.81466, 0.793893, 0.77232, 0.75, 0.726995, 0.703368, 0.679184, 0.654508, 0.62941, 0.603956, 0.578217, 0.552264, 0.526168, 0.5, 0.473832, 0.447736, 0.421783, 0.396044, 0.37059, 0.345492, 0.320816, 0.296632, 0.273005, 0.25, 0.22768, 0.206107, 0.18534, 0.165435, 0.146447, 0.128428, 0.111427, 0.095492, 0.080665, 0.066987, 0.054497, 0.043227, 0.03321, 0.024472, 0.017037, 0.010926, 0.006156, 0.002739, 0.000685, 0.0, 0.000685, 0.002739, 0.006156, 0.010926, 0.017037, 0.024472, 0.03321, 0.043227, 0.054497, 0.066987, 0.080665, 0.095492, 0.111427, 0.128428, 0.146447, 0.165435, 0.18534, 0.206107, 0.22768, 0.25, 0.273005, 0.296632, 0.320816, 0.345492, 0.37059, 0.396044, 0.421783, 0.447736, 0.473832, 0.5, 0.526168, 0.552264, 0.578217, 0.603956, 0.62941, 0.654508, 0.679184, 0.703368, 0.726995, 0.75, 0.77232, 0.793893, 0.81466, 0.834565, 0.853553, 0.871572, 0.888573, 0.904508, 0.919335, 0.933013, 0.945503, 0.956773, 0.96679, 0.975528, 0.982963, 0.989074, 0.993844, 0.997261};
y={-0.0, 0.000397, 0.000891, 0.001577, 0.002449, 0.003501, 0.004725, 0.006112, 0.007651, 0.009332, 0.011142, 0.013071, 0.015105, 0.017232, 0.019438, 0.021712, 0.024038, 0.026405, 0.028798, 0.031204, 0.03361, 0.036002, 0.038366, 0.040686, 0.042949, 0.04514, 0.047242, 0.049241, 0.051119, 0.052862, 0.054451, 0.055872, 0.057108, 0.058144, 0.058965, 0.059557, 0.059907, 0.060005, 0.059841, 0.059407, 0.058699, 0.057712, 0.056446, 0.054901, 0.053083, 0.050995, 0.048648, 0.046049, 0.04321, 0.040145, 0.036867, 0.033389, 0.029726, 0.025893, 0.021904, 0.01777, 0.013503, 0.009114, 0.004611, 0.0, -0.004611, -0.009114, -0.013503, -0.01777, -0.021904, -0.025893, -0.029726, -0.033389, -0.036867, -0.040145, -0.04321, -0.046049, -0.048648, -0.050995, -0.053083, -0.054901, -0.056446, -0.057712, -0.058699, -0.059407, -0.059841, -0.060005, -0.059907, -0.059557, -0.058965, -0.058144, -0.057108, -0.055872, -0.054451, -0.052862, -0.051119, -0.049241, -0.047242, -0.04514, -0.042949, -0.040686, -0.038366, -0.036002, -0.03361, -0.031204, -0.028798, -0.026405, -0.024038, -0.021712, -0.019438, -0.017232, -0.015105, -0.013071, -0.011142, -0.009332, -0.007651, -0.006112, -0.004725, -0.003501, -0.002449, -0.001577, -0.000891, -0.000397};
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
