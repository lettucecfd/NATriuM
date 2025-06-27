//FILE newNonUniCoarsen2
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
n_coarsen    = 2;
n_points_orig = 197;
n_points     = 99;
point_id_top = 34;  // 68
point_id_bot = 116;  // (n_points_orig+1) - (point_id_top-2) - 1;
point_id_front = 100;
n_around     = 50; // 50/n_coarsen for plot
n_foil       = 2;
n_inlet      = 51 - point_id_top;
channel_l    = 1.5;
channel_h    = inlet_r;
n_channel    = 200;  // 300/n_coarsen for plot
n_outlet_center = n_channel - point_id_top + 1;
progression_around = 1.05;
progression_sponge_front = 1.11;
progression_sponge_back = 1.1;
n_sponge_front= 29;
n_sponge_back= 26;


/// FOIL
x = {1.0, 0.999748, 0.998993, 0.997736, 0.995977, 0.993719, 0.990964, 0.987715, 0.983974, 0.979746, 0.975036, 0.969846, 0.964184, 0.958054, 0.951463, 0.944418, 0.936925, 0.928992, 0.920627, 0.911838, 0.902635, 0.893027, 0.883022, 0.872632, 0.861867, 0.850737, 0.839255, 0.82743, 0.815276, 0.802805, 0.790028, 0.77696, 0.763613, 0.75, 0.736136, 0.722033, 0.707708, 0.693173, 0.678443, 0.663534, 0.64846, 0.633237, 0.617879, 0.602403, 0.586824, 0.571157, 0.555419, 0.539625, 0.523791, 0.507933, 0.492067, 0.476209, 0.460375, 0.444581, 0.428843, 0.413176, 0.397597, 0.382121, 0.366763, 0.35154, 0.336466, 0.321557, 0.292292, 0.277967, 0.263864, 0.25, 0.236387, 0.22304, 0.209972, 0.197195, 0.184724, 0.17257, 0.160745, 0.149263, 0.138133, 0.127368, 0.116978, 0.106973, 0.0973649, 0.0881617, 0.0793732, 0.0710083, 0.0630753, 0.0555823, 0.0485367, 0.0419458, 0.035816, 0.0301537, 0.0249644, 0.0202535, 0.0160256, 0.0122851, 0.00903565, 0.00628056, 0.00402259, 0.00226404, 0.00100666, 0.000251729, 0.0, 0.000251729, 0.00100666, 0.00226404, 0.00402259, 0.00628056, 0.00903565, 0.0122851, 0.0160256, 0.0202535, 0.0249644, 0.0301537, 0.035816, 0.0419458, 0.0485367, 0.0555823, 0.0630753, 0.0710083, 0.0793732, 0.0881617, 0.0973649, 0.106973, 0.116978, 0.127368, 0.138133, 0.149263, 0.160745, 0.17257, 0.184724, 0.197195, 0.209972, 0.22304, 0.236387, 0.25, 0.263864, 0.277967, 0.306827, 0.321557, 0.336466, 0.35154, 0.366763, 0.382121, 0.397597, 0.413176, 0.428843, 0.444581, 0.460375, 0.476209, 0.492067, 0.507933, 0.523791, 0.539625, 0.555419, 0.571157, 0.586824, 0.602403, 0.617879, 0.633237, 0.64846, 0.663534, 0.678443, 0.693173, 0.707708, 0.722033, 0.736136, 0.75, 0.763613, 0.77696, 0.790028, 0.802805, 0.815276, 0.82743, 0.839255, 0.850737, 0.861867, 0.872632, 0.883022, 0.893027, 0.902635, 0.911838, 0.920627, 0.928992, 0.936925, 0.944418, 0.951463, 0.958054, 0.964184, 0.969846, 0.975036, 0.979746, 0.983974, 0.987715, 0.990964, 0.993719, 0.995977, 0.997736, 0.998993, 0.999748};
y = {-1.66533e-17, 3.65828e-05, 0.000146223, 0.000328595, 0.00058316, 0.00090917, 0.00130567, 0.00177151, 0.00230534, 0.00290565, 0.00357073, 0.00429874, 0.00508767, 0.00593537, 0.00683958, 0.00779793, 0.00880793, 0.00986703, 0.0109726, 0.0121219, 0.0133121, 0.0145406, 0.0158044, 0.0171006, 0.0184265, 0.0197789, 0.0211551, 0.0225521, 0.0239669, 0.0253966, 0.0268382, 0.0282887, 0.0297451, 0.0312044, 0.0326635, 0.0341194, 0.035569, 0.0370091, 0.0384366, 0.0398482, 0.0412407, 0.0426107, 0.0439549, 0.0452698, 0.0465521, 0.0477982, 0.0490045, 0.0501675, 0.0512835, 0.0523489, 0.0533601, 0.0543134, 0.0552053, 0.0560321, 0.0567903, 0.0574765, 0.0580873, 0.0586193, 0.0590696, 0.059435, 0.0597128, 0.0599003, 0.0599951, 0.0598982, 0.0597028, 0.0594075, 0.0590111, 0.0585128, 0.0579121, 0.0572087, 0.0564028, 0.0554948, 0.0544854, 0.0533756, 0.0521668, 0.0508605, 0.0494586, 0.0479633, 0.0463768, 0.0447018, 0.042941, 0.0410972, 0.0391735, 0.037173, 0.0350989, 0.0329544, 0.0307426, 0.0284669, 0.0261302, 0.0237357, 0.0212862, 0.0187844, 0.0162331, 0.0136345, 0.0109908, 0.008304, 0.0055757, 0.00280732, 0.0, -0.00280732, -0.0055757, -0.008304, -0.0109908, -0.0136345, -0.0162331, -0.0187844, -0.0212862, -0.0237357, -0.0261302, -0.0284669, -0.0307426, -0.0329544, -0.0350989, -0.037173, -0.0391735, -0.0410972, -0.042941, -0.0447018, -0.0463768, -0.0479633, -0.0494586, -0.0508605, -0.0521668, -0.0533756, -0.0544854, -0.0554948, -0.0564028, -0.0572087, -0.0579121, -0.0585128, -0.0590111, -0.0594075, -0.0597028, -0.0598982, -0.0599951, -0.0599003, -0.0597128, -0.059435, -0.0590696, -0.0586193, -0.0580873, -0.0574765, -0.0567903, -0.0560321, -0.0552053, -0.0543134, -0.0533601, -0.0523489, -0.0512835, -0.0501675, -0.0490045, -0.0477982, -0.0465521, -0.0452698, -0.0439549, -0.0426107, -0.0412407, -0.0398482, -0.0384366, -0.0370091, -0.035569, -0.0341194, -0.0326635, -0.0312044, -0.0297451, -0.0282887, -0.0268382, -0.0253966, -0.0239669, -0.0225521, -0.0211551, -0.0197789, -0.0184265, -0.0171006, -0.0158044, -0.0145406, -0.0133121, -0.0121219, -0.0109726, -0.00986703, -0.00880793, -0.00779793, -0.00683958, -0.00593537, -0.00508767, -0.00429874, -0.00357073, -0.00290565, -0.00230534, -0.00177151, -0.00130567, -0.00090917, -0.00058316, -0.000328595, -0.000146223, -3.65828e-05};

nx = n_points;//#x[]/n_coarsen-1;


For i In {1:nx/2}
  Point(i) = {x[n_coarsen*(i-1)], y[n_coarsen*(i-1)], 0, 1};
EndFor
Point(point_id_front) =  {x[98], y[98], 0, 1};
For i In {101:99+nx/2}
  Point(i) = {x[n_coarsen*(i-51)], y[n_coarsen*(i-51)], 0, 1};
EndFor



For i In {1:48}
  Line(i) = {i, i+1};
EndFor
Line(49) = {49, point_id_front};
For i In {50:nx-2}
  iPoint = 100+ i-50;
  Line(i) = {iPoint, iPoint+1};
EndFor
Line(nx-1) = {148, 1};



Transfinite Curve {1:nx-1} = n_foil;


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

Curve Loop (200) = {-203, 201, point_id_top:49, 200}; // inlet loop top
Curve Loop (201) = {-200, 50:65, -202, -204}; // inlet loop bottom
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
Curve Loop(223) = {202, 66:98, 210, -212, -221}; // channel bottom
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
Physical Curve(300) = {233, 204, 203, 232};  // "Inlet", 
// TODO: only inlet in channel, not sponge Physical Curve(300) = {204, 203};  // "Inlet", 
//Physical Curve("Sponge", 301) = {234, 235};
Physical Curve(302) = {231, 212, 211, 230};  // "Outlet", 
Physical Curve(303) = {1:nx-1};  // "BB_BC", 
Physical Surface(236) = {1, 2, 6, 5, 7, 8};

Mesh 2;
//RecombineMesh;
Save "NACA0012_0deg.msh";