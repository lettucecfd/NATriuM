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
n_coarsen    = 1;
n_points     = 286;  // 197/n_coarsen;
point_id_top = 105;  // 76/n_coarsen;
point_id_bot = 185;  // (n_points+1) - (point_id_top-2);
point_id_front = 145;  // 99/n_coarsen + 1;
n_around     = 71;  // 100/n_coarsen; // 50/n_coarsen for plot
n_foil       = 2;
n_inlet      = point_id_front - point_id_top + 1;  // (n_coarsen+1)/2;  // +1 @1; +2 @3
channel_l    = 1.5;
channel_h    = inlet_r;
n_channel    = 80;  // 400/n_coarsen;  // 300/n_coarsen for plot
n_outlet_center = 901 - point_id_top + 1;  // (n_coarsen+1)/2;  // +1 @1; +2 @3
progression_around = 1.05;
progression_sponge = 1.05;
n_sponge_front= 30;
n_sponge_back= 25;


/// FOIL
x={1.0, 1.0, 0.99, 0.99, 0.98, 0.98, 0.97, 0.97, 0.96, 0.96, 0.95, 0.95, 0.94, 0.94, 0.93, 0.93, 0.92, 0.92, 0.91, 0.91, 0.9, 0.9, 0.89, 0.89, 0.88, 0.88, 0.87, 0.87, 0.86, 0.86, 0.85, 0.85, 0.84, 0.84, 0.83, 0.83, 0.82, 0.82, 0.81, 0.81, 0.8, 0.8, 0.79, 0.79, 0.78, 0.78, 0.77, 0.77, 0.76, 0.76, 0.75, 0.75, 0.74, 0.74, 0.73, 0.73, 0.72, 0.72, 0.71, 0.7, 0.69, 0.68, 0.67, 0.66, 0.65, 0.64, 0.63, 0.62, 0.61, 0.6, 0.59, 0.58, 0.57, 0.56, 0.55, 0.54, 0.53, 0.52, 0.51, 0.5, 0.49, 0.48, 0.47, 0.46, 0.45, 0.44, 0.43, 0.42, 0.41, 0.4, 0.39, 0.38, 0.37, 0.36, 0.35, 0.34, 0.33, 0.32, 0.31, 0.3, 0.29, 0.28, 0.27, 0.26, 0.25, 0.24, 0.23, 0.22, 0.21, 0.2, 0.19, 0.18, 0.17, 0.16, 0.15, 0.15, 0.14, 0.14, 0.13, 0.13, 0.12, 0.12, 0.11, 0.11, 0.1, 0.1, 0.09, 0.09, 0.08, 0.08, 0.07, 0.07, 0.06, 0.06, 0.05, 0.05, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.01, 0.01, 0.0, 0.0, 0.01, 0.01, 0.02, 0.02, 0.03, 0.03, 0.04, 0.04, 0.05, 0.05, 0.06, 0.06, 0.07, 0.07, 0.08, 0.08, 0.09, 0.09, 0.1, 0.1, 0.11, 0.11, 0.12, 0.12, 0.13, 0.13, 0.14, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.4, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.5, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.6, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.7, 0.71, 0.71, 0.72, 0.72, 0.73, 0.73, 0.74, 0.74, 0.75, 0.75, 0.76, 0.76, 0.77, 0.77, 0.78, 0.78, 0.79, 0.79, 0.8, 0.8, 0.81, 0.81, 0.82, 0.82, 0.83, 0.83, 0.84, 0.84, 0.85, 0.85, 0.86, 0.86, 0.87, 0.87, 0.88, 0.88, 0.89, 0.89, 0.9, 0.9, 0.91, 0.91, 0.92, 0.92, 0.93, 0.93, 0.94, 0.94, 0.95, 0.95, 0.96, 0.96, 0.97, 0.97, 0.98, 0.98, 0.99};
y={-0.0, 0.001444, 0.001444, 0.00287, 0.00287, 0.004277, 0.004277, 0.005667, 0.005667, 0.007039, 0.007039, 0.008395, 0.008395, 0.009733, 0.009733, 0.011055, 0.011055, 0.012361, 0.012361, 0.01365, 0.01365, 0.014925, 0.014925, 0.016183, 0.016183, 0.017426, 0.017426, 0.018655, 0.018655, 0.019868, 0.019868, 0.021066, 0.021066, 0.02225, 0.02225, 0.02342, 0.02342, 0.024575, 0.024575, 0.025715, 0.025715, 0.026841, 0.026841, 0.027953, 0.027953, 0.029051, 0.029051, 0.030135, 0.030135, 0.031204, 0.031204, 0.03226, 0.03226, 0.0333, 0.0333, 0.034327, 0.034327, 0.035339, 0.035339, 0.036337, 0.037319, 0.038287, 0.03924, 0.040178, 0.0411, 0.042007, 0.042897, 0.043772, 0.044629, 0.04547, 0.046294, 0.0471, 0.047888, 0.048658, 0.049409, 0.05014, 0.050852, 0.051543, 0.052213, 0.052862, 0.053488, 0.054091, 0.05467, 0.055226, 0.055756, 0.05626, 0.056737, 0.057186, 0.057607, 0.057998, 0.058358, 0.058686, 0.058981, 0.059242, 0.059467, 0.059655, 0.059805, 0.059915, 0.059983, 0.060007, 0.059986, 0.059918, 0.0598, 0.059631, 0.059408, 0.059127, 0.058787, 0.058383, 0.057914, 0.057373, 0.056759, 0.056065, 0.055286, 0.054417, 0.053451, 0.052379, 0.052379, 0.051193, 0.051193, 0.049882, 0.049882, 0.048432, 0.048432, 0.046828, 0.046828, 0.045049, 0.045049, 0.043072, 0.043072, 0.040863, 0.040863, 0.038376, 0.038376, 0.035547, 0.035547, 0.032277, 0.032277, 0.028401, 0.028401, 0.023598, 0.023598, 0.017037, 0.017037, 0.0, 0.0, -0.017037, -0.017037, -0.023598, -0.023598, -0.028401, -0.028401, -0.032277, -0.032277, -0.035547, -0.035547, -0.038376, -0.038376, -0.040863, -0.040863, -0.043072, -0.043072, -0.045049, -0.045049, -0.046828, -0.046828, -0.048432, -0.048432, -0.049882, -0.049882, -0.051193, -0.051193, -0.052379, -0.052379, -0.053451, -0.053451, -0.054417, -0.055286, -0.056065, -0.056759, -0.057373, -0.057914, -0.058383, -0.058787, -0.059127, -0.059408, -0.059631, -0.0598, -0.059918, -0.059986, -0.060007, -0.059983, -0.059915, -0.059805, -0.059655, -0.059467, -0.059242, -0.058981, -0.058686, -0.058358, -0.057998, -0.057607, -0.057186, -0.056737, -0.05626, -0.055756, -0.055226, -0.05467, -0.054091, -0.053488, -0.052862, -0.052213, -0.051543, -0.050852, -0.05014, -0.049409, -0.048658, -0.047888, -0.0471, -0.046294, -0.04547, -0.044629, -0.043772, -0.042897, -0.042007, -0.0411, -0.040178, -0.03924, -0.038287, -0.037319, -0.036337, -0.035339, -0.034327, -0.034327, -0.0333, -0.0333, -0.03226, -0.03226, -0.031204, -0.031204, -0.030135, -0.030135, -0.029051, -0.029051, -0.027953, -0.027953, -0.026841, -0.026841, -0.025715, -0.025715, -0.024575, -0.024575, -0.02342, -0.02342, -0.02225, -0.02225, -0.021066, -0.021066, -0.019868, -0.019868, -0.018655, -0.018655, -0.017426, -0.017426, -0.016183, -0.016183, -0.014925, -0.014925, -0.01365, -0.01365, -0.012361, -0.012361, -0.011055, -0.011055, -0.009733, -0.009733, -0.008395, -0.008395, -0.007039, -0.007039, -0.005667, -0.005667, -0.004277, -0.004277, -0.00287, -0.00287, -0.001444, -0.001444};
nx = n_points;//#x[]/n_coarsen-1;

For i In {0: nx}
  Point(i+1) = {x[n_coarsen*i], y[n_coarsen*i], 0, 0.004};
EndFor

For i In {1: nx}
  Line(i) = {i, i+1};
EndFor
Line(nx+1) = {nx+1, 1};

For i In {1:57:2}
  Transfinite Curve {i} = 2 Using Progression 1;
EndFor
For i In {115:137:2}
  Transfinite Curve {i} = 2 Using Progression 1;
EndFor
For i In {153:173:2}
  Transfinite Curve {i} = 2 Using Progression 1;
EndFor
For i In {231:285:2}
  Transfinite Curve {i} = 2 Using Progression 1;
EndFor


/// INLET
Point(400) = {0, inlet_r, 0, size_in_out};       // inlet top
Point(401) = {0, -inlet_r, 0, size_in_out};      // inlet bottom
Point(402) = {inlet_c-inlet_front, 0, 0, size_in_out};     // inlet front
Point(403) = {inlet_c, 0, 0};                          // inlet center
Line(500)  = {point_id_front, 402};                               // inlet front line
Line(501)  = {400, point_id_top};                      // inlet top line
Line(502)  = {401, point_id_bot};                      // inlet bottom line
Ellipse(503) = {400, 403, 400, 402};                   // inlet circle line top
Ellipse(504) = {401, 403, 401, 402};                   // inlet circle line bottom
Transfinite Curve {500, -501} = n_around Using Progression progression_around;  // inlet lines points
Transfinite Curve {-502} = n_around-1 Using Progression progression_around;  // inlet lines points
Transfinite Curve {503, 504} = n_inlet Using Progression 1;        // inlet circle points

Curve Loop (600) = {-503, 501, point_id_top:(point_id_front-1), 500};   // inlet loop top
Curve Loop (601) = {-500, point_id_front:(point_id_bot-1), -502, 504};  // inlet loop bottom
Plane Surface(1) = {600};                                // inlet surface top
Plane Surface(2) = {601};                                // inlet surface bottom
//Transfinite Surface {1} = {402, point_id_front, point_id_top, 400};           // inlet surface top
//Transfinite Surface {2} = {402, 401, point_id_bot, point_id_front};          // inlet surface bottom


/// OUTLET
Point(410) = {outlet_c, outlet_h, 0, size_in_out};     // outlet top
Point(411) = {outlet_c, -outlet_h, 0, size_in_out};    // outlet bottom
Point(412) = {outlet_c, 0, 0, size_foil};              // outlet center
Line(510)  = {1, 412};                                 // outlet center line
Line(511)  = {410, 412};                               // outlet end top line
Line(512)  = {411, 412};                               // outlet end bottom line
Transfinite Curve {510} = n_outlet_center Using Progression 1;                // outlet center lines points
Transfinite Curve {-511} = n_around+0 Using Progression progression_around;   // outlet end lines points top
Transfinite Curve {-512} = n_around+0 Using Progression progression_around;   // outlet end lines points bottom


/// CHANNEL
Line(520)  = {400, 410};                                  // channel top line
Line(521)  = {401, 411};                                  // channel bottom line
Transfinite Curve {520} = n_channel Using Progression 1;  // channel top points
Transfinite Curve {521} = n_channel Using Progression 1;  // channel bottom points


/// OUTLET / CHANNEL SURFACES
Curve Loop(622) = {-501, 520, 511, -510, 1:(point_id_top-1)};   // channel top
Curve Loop(623) = {502, point_id_bot:(nx+1), 510, -512, -521};  // channel bottom
Plane Surface(5) = {622};                                       // channel top
Plane Surface(6) = {623};                                       // channel bottom
//Transfinite Surface {5} = {400, point_id_top, 412, 410};      // channel surface top
//Transfinite Surface {6} = {point_id_bot, 401, 411, 412};      // channel surface bottom


/// SPONGE
Point(430) = {outlet_c, sponge_h, 0, size_sponge};      // sponge top back
Point(431) = {outlet_c, -sponge_h, 0, size_sponge};     // sponge bottom back
Line(530)  = {410, 430};                                // sponge outlet top line
Line(531)  = {411, 431};                                // sponge outlet bottom line
Line(532)  = {400, 430};                                // sponge top line
Line(533)  = {401, 431};                                // sponge bottom line
Curve Loop(630) = {532, -530, -520};                    // sponge surface top
Curve Loop(631) = {521, 531, -533};                     // sponge surface bottom
Transfinite Curve {532, 533} = n_sponge_front Using Progression progression_sponge;   // sponge front/diagonal lines
Transfinite Curve {530, 531} = n_sponge_back Using Progression progression_sponge;    // sponge back lines
Plane Surface(7) = {630};                               // sponge top
Plane Surface(8) = {631};                               // sponge bottom

/// MESH SIZES
Mesh.Algorithm = 6;
Mesh.RecombineAll = 1;
Mesh.RecombinationAlgorithm = 1; // or 3; to leave no triangles

/// BOUNDARIES
Physical Curve("Inlet", 300) = {533, 504, 503, 532};
//Physical Curve("Sponge", 301) = {534, 535};
Physical Curve("Outlet", 302) = {531, 512, 511, 530};
Physical Curve("BB_BC", 303) = {1:(nx+1)};
Physical Surface(236) = {1, 2, 6, 5, 7, 8};

Mesh 2;
Save "NACA0012_0deg.msh";

//Point(432) = {-2.7, 0.1, 0, 1.0};
//Point(433) = {-0.1, 0.6, 0, 1.0};
//Point(434) = {-2.7, -0.1, 0, 1.0};
//Line(534) = {432, 115};
//Line(535) = {174, 434};
//Point(435) = {6, -0.05, 0, 1.0};
//Point(436) = {6, 0.05, 0, 1.0};
//Line(536) = {436, 58};
//Line(537) = {231, 435};
