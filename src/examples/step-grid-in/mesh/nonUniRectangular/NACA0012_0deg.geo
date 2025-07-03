//FILE newNonUni 0deg res100

inlet_r      = 4;
inlet_front  = 3;
inlet_c      = 0;
outlet_c     = 6;
outlet_h     = inlet_r;
sponge_h     = inlet_r * 3;
size_foil    = 0.1;
size_in_out  = 1;
size_sponge  = 1;
point_id_top = 64;
point_id_front = 100;
point_id_bot = 136;
n_around     = 70.0;
n_foil       = 2;  // number of points per foil section; 2 for just the section, 3 to split once, ...
n_inlet      = point_id_front - point_id_top + 1;  // top and bottom, each; 31 default points; additional 31-1 (30 segments) points for each n_foil above 2
channel_l    = 1.5;
channel_h    = inlet_r;
n_channel    = 500;  // number of points on top and bottom
n_outlet_center = n_channel - point_id_top + 1;
progression_around = 1.05;

x={1.0, 0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.91, 0.9, 0.89, 0.88, 0.87, 0.86, 0.85, 0.84, 0.83, 0.82, 0.81, 0.8, 0.79, 0.78, 0.77, 0.76, 0.75, 0.74, 0.73, 0.72, 0.71, 0.7, 0.69, 0.68, 0.67, 0.66, 0.65, 0.64, 0.63, 0.62, 0.61, 0.6, 0.59, 0.58, 0.57, 0.56, 0.55, 0.54, 0.53, 0.515705, 0.5, 0.484295, 0.468605, 0.452946, 0.437333, 0.421783, 0.406309, 0.390928, 0.375655, 0.360504, 0.345492, 0.330631, 0.315938, 0.301426, 0.28711, 0.273005, 0.259123, 0.245479, 0.232087, 0.218958, 0.206107, 0.193546, 0.181288, 0.169344, 0.157726, 0.146447, 0.135516, 0.124944, 0.114743, 0.104922, 0.095492, 0.08646, 0.077836, 0.069629, 0.061847, 0.054497, 0.047586, 0.041123, 0.035112, 0.02956, 0.024472, 0.019853, 0.015708, 0.012042, 0.008856, 0.006156, 0.003943, 0.002219, 0.000987, 0.000247, 0.0, 0.000247, 0.000987, 0.002219, 0.003943, 0.006156, 0.008856, 0.012042, 0.015708, 0.019853, 0.024472, 0.02956, 0.035112, 0.041123, 0.047586, 0.054497, 0.061847, 0.069629, 0.077836, 0.08646, 0.095492, 0.104922, 0.114743, 0.124944, 0.135516, 0.146447, 0.157726, 0.169344, 0.181288, 0.193546, 0.206107, 0.218958, 0.232087, 0.245479, 0.259123, 0.273005, 0.28711, 0.301426, 0.315938, 0.330631, 0.345492, 0.360504, 0.375655, 0.390928, 0.406309, 0.421783, 0.437333, 0.452946, 0.468605, 0.484295, 0.5, 0.515705, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.8, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.9, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99};
y={-0.0, 0.001444, 0.00287, 0.004277, 0.005667, 0.007039, 0.008395, 0.009733, 0.011055, 0.012361, 0.01365, 0.014925, 0.016183, 0.017426, 0.018655, 0.019868, 0.021066, 0.02225, 0.02342, 0.024575, 0.025715, 0.026841, 0.027953, 0.029051, 0.030135, 0.031204, 0.03226, 0.0333, 0.034327, 0.035339, 0.036337, 0.037319, 0.038287, 0.03924, 0.040178, 0.0411, 0.042007, 0.042897, 0.043772, 0.044629, 0.04547, 0.046294, 0.0471, 0.047888, 0.048658, 0.049409, 0.05014, 0.050852, 0.051833, 0.052862, 0.053835, 0.054749, 0.055602, 0.05639, 0.057108, 0.057755, 0.058326, 0.058819, 0.05923, 0.059557, 0.059797, 0.059947, 0.060006, 0.059971, 0.059841, 0.059614, 0.059288, 0.058863, 0.058338, 0.057712, 0.056986, 0.056159, 0.055232, 0.054206, 0.053083, 0.051862, 0.050546, 0.049138, 0.047638, 0.046049, 0.044374, 0.042615, 0.040776, 0.038859, 0.036867, 0.034803, 0.032671, 0.030473, 0.028213, 0.025893, 0.023517, 0.021088, 0.018607, 0.016078, 0.013503, 0.010884, 0.008223, 0.005521, 0.002779, 0.0, -0.002779, -0.005521, -0.008223, -0.010884, -0.013503, -0.016078, -0.018607, -0.021088, -0.023517, -0.025893, -0.028213, -0.030473, -0.032671, -0.034803, -0.036867, -0.038859, -0.040776, -0.042615, -0.044374, -0.046049, -0.047638, -0.049138, -0.050546, -0.051862, -0.053083, -0.054206, -0.055232, -0.056159, -0.056986, -0.057712, -0.058338, -0.058863, -0.059288, -0.059614, -0.059841, -0.059971, -0.060006, -0.059947, -0.059797, -0.059557, -0.05923, -0.058819, -0.058326, -0.057755, -0.057108, -0.05639, -0.055602, -0.054749, -0.053835, -0.052862, -0.051833, -0.050852, -0.05014, -0.049409, -0.048658, -0.047888, -0.0471, -0.046294, -0.04547, -0.044629, -0.043772, -0.042897, -0.042007, -0.0411, -0.040178, -0.03924, -0.038287, -0.037319, -0.036337, -0.035339, -0.034327, -0.0333, -0.03226, -0.031204, -0.030135, -0.029051, -0.027953, -0.026841, -0.025715, -0.024575, -0.02342, -0.02225, -0.021066, -0.019868, -0.018655, -0.017426, -0.016183, -0.014925, -0.01365, -0.012361, -0.011055, -0.009733, -0.008395, -0.007039, -0.005667, -0.004277, -0.00287, -0.001444};
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

/// MESH SIZES
//Mesh.ElementOrder = 1;
//Mesh.Algorithm = 6;
Mesh.RecombineAll = 1;
Mesh.SubdivisionAlgorithm = -1;
Mesh.RecombinationAlgorithm = 1; // or 3; to leave no triangles

/// BOUNDARIES
Physical Curve(300) = {220, 204, 203, 221};  // "Inlet", 
Physical Curve(302) = {212, 211};  // "Outlet", 
Physical Curve(303) = {1:nx};  // "BB_BC", 
Physical Surface(236) = {1, 2, 6, 5};

Mesh 2;
Save "NACA0012_0deg.msh";