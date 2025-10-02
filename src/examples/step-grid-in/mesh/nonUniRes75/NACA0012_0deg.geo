//FILE newNonUni 0deg res75

inlet_r      = 4;
inlet_front  = 3;
inlet_c      = 0.3;
outlet_c     = 6;
outlet_h     = inlet_r;
sponge_h     = inlet_r * 3;
size_foil    = 0.1;
size_in_out  = 1;
size_sponge  = 1;
point_id_top = 49;
point_id_front = 75;
point_id_bot = 101;
n_around     = 75;
n_foil       = 2;  // number of points per foil section; 2 for just the section, 3 to split once, ...
n_inlet      = point_id_front - point_id_top + 1;  // top and bottom, each; 31 default points; additional 31-1 (30 segments) points for each n_foil above 2
channel_l    = 1.5;
channel_h    = inlet_r;
n_channel    = 262;  // number of points on top and bottom
n_outlet_center = n_channel - point_id_top + 1;
progression_around = 1.05;
progression_sponge_front = 1.05;
progression_sponge_back = 1.05;
n_sponge_front= 57;
n_sponge_back= 50;

x={1.0, 0.998246, 0.996057, 0.992998, 0.989074, 0.984292, 0.97866, 0.972188, 0.964888, 0.956773, 0.947856, 0.938153, 0.927682, 0.916461, 0.904508, 0.891847, 0.878498, 0.864484, 0.849832, 0.834565, 0.818712, 0.8023, 0.785357, 0.767913, 0.75, 0.731648, 0.71289, 0.693758, 0.674286, 0.654508, 0.63446, 0.614175, 0.593691, 0.573042, 0.552264, 0.531395, 0.510471, 0.489529, 0.468605, 0.447736, 0.426958, 0.406309, 0.385825, 0.36554, 0.345492, 0.325714, 0.306242, 0.28711, 0.268352, 0.25, 0.232087, 0.214643, 0.1977, 0.181288, 0.165435, 0.150168, 0.135516, 0.121502, 0.108153, 0.095492, 0.083539, 0.072318, 0.061847, 0.052144, 0.043227, 0.035112, 0.027812, 0.02134, 0.015708, 0.010926, 0.007002, 0.003943, 0.001754, 0.000439, 0.0, 0.000439, 0.001754, 0.003943, 0.007002, 0.010926, 0.015708, 0.02134, 0.027812, 0.035112, 0.043227, 0.052144, 0.061847, 0.072318, 0.083539, 0.095492, 0.108153, 0.121502, 0.135516, 0.150168, 0.165435, 0.181288, 0.1977, 0.214643, 0.232087, 0.25, 0.268352, 0.28711, 0.306242, 0.325714, 0.345492, 0.36554, 0.385825, 0.406309, 0.426958, 0.447736, 0.468605, 0.489529, 0.510471, 0.531395, 0.552264, 0.573042, 0.593691, 0.614175, 0.63446, 0.654508, 0.674286, 0.693758, 0.71289, 0.731648, 0.75, 0.767913, 0.785357, 0.8023, 0.818712, 0.834565, 0.849832, 0.864484, 0.878498, 0.891847, 0.904508, 0.916461, 0.927682, 0.938153, 0.947856, 0.956773, 0.964888, 0.972188, 0.97866, 0.984292, 0.989074, 0.992998, 0.996057, 0.998246};
y={-0.0, 0.000255, 0.000572, 0.001013, 0.001577, 0.00226, 0.003059, 0.003971, 0.00499, 0.006112, 0.007331, 0.008643, 0.010041, 0.011519, 0.013071, 0.01469, 0.016371, 0.018106, 0.019888, 0.021712, 0.023569, 0.025454, 0.02736, 0.029279, 0.031204, 0.03313, 0.035048, 0.036952, 0.038834, 0.040686, 0.042502, 0.044273, 0.045992, 0.047651, 0.049241, 0.050754, 0.052182, 0.053517, 0.054749, 0.055872, 0.056877, 0.057755, 0.058499, 0.059102, 0.059557, 0.059857, 0.059997, 0.059971, 0.059776, 0.059407, 0.058863, 0.05814, 0.057239, 0.056159, 0.054901, 0.053468, 0.051862, 0.050087, 0.048148, 0.046049, 0.043797, 0.041398, 0.038859, 0.036186, 0.033389, 0.030473, 0.027446, 0.024315, 0.021088, 0.01777, 0.014367, 0.010884, 0.007327, 0.003697, 0.0, -0.003697, -0.007327, -0.010884, -0.014367, -0.01777, -0.021088, -0.024315, -0.027446, -0.030473, -0.033389, -0.036186, -0.038859, -0.041398, -0.043797, -0.046049, -0.048148, -0.050087, -0.051862, -0.053468, -0.054901, -0.056159, -0.057239, -0.05814, -0.058863, -0.059407, -0.059776, -0.059971, -0.059997, -0.059857, -0.059557, -0.059102, -0.058499, -0.057755, -0.056877, -0.055872, -0.054749, -0.053517, -0.052182, -0.050754, -0.049241, -0.047651, -0.045992, -0.044273, -0.042502, -0.040686, -0.038834, -0.036952, -0.035048, -0.03313, -0.031204, -0.029279, -0.02736, -0.025454, -0.023569, -0.021712, -0.019888, -0.018106, -0.016371, -0.01469, -0.013071, -0.011519, -0.010041, -0.008643, -0.007331, -0.006112, -0.00499, -0.003971, -0.003059, -0.00226, -0.001577, -0.001013, -0.000572, -0.000255};
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
