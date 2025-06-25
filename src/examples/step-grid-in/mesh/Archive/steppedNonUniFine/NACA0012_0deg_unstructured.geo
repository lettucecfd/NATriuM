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

x = {1.0, 0.999748, 0.998993, 0.997736, 0.995977, 0.993719, 0.990964, 0.987715, 0.983974, 0.979746, 0.975036, 0.969846, 0.964184, 0.958054, 0.951463, 0.944418, 0.936925, 0.928992, 0.920627, 0.911838, 0.902635, 0.893027, 0.883022, 0.872632, 0.861867, 0.850737, 0.839255, 0.82743, 0.815276, 0.802805, 0.790028, 0.77696, 0.763613, 0.75, 0.736136, 0.722033, 0.707708, 0.693173, 0.678443, 0.663534, 0.64846, 0.633237, 0.617879, 0.602403, 0.586824, 0.571157, 0.555419, 0.539625, 0.523791, 0.507933, 0.492067, 0.476209, 0.460375, 0.444581, 0.428843, 0.413176, 0.397597, 0.382121, 0.366763, 0.35154, 0.336466, 0.321557, 0.292292, 0.277967, 0.263864, 0.25, 0.236387, 0.22304, 0.209972, 0.197195, 0.184724, 0.17257, 0.160745, 0.149263, 0.138133, 0.127368, 0.116978, 0.106973, 0.0973649, 0.0881617, 0.0793732, 0.0710083, 0.0630753, 0.0555823, 0.0485367, 0.0419458, 0.035816, 0.0301537, 0.0249644, 0.0202535, 0.0160256, 0.0122851, 0.00903565, 0.00628056, 0.00402259, 0.00226404, 0.00100666, 0.000251729, 0.0, 0.000251729, 0.00100666, 0.00226404, 0.00402259, 0.00628056, 0.00903565, 0.0122851, 0.0160256, 0.0202535, 0.0249644, 0.0301537, 0.035816, 0.0419458, 0.0485367, 0.0555823, 0.0630753, 0.0710083, 0.0793732, 0.0881617, 0.0973649, 0.106973, 0.116978, 0.127368, 0.138133, 0.149263, 0.160745, 0.17257, 0.184724, 0.197195, 0.209972, 0.22304, 0.236387, 0.25, 0.263864, 0.277967, 0.306827, 0.321557, 0.336466, 0.35154, 0.366763, 0.382121, 0.397597, 0.413176, 0.428843, 0.444581, 0.460375, 0.476209, 0.492067, 0.507933, 0.523791, 0.539625, 0.555419, 0.571157, 0.586824, 0.602403, 0.617879, 0.633237, 0.64846, 0.663534, 0.678443, 0.693173, 0.707708, 0.722033, 0.736136, 0.75, 0.763613, 0.77696, 0.790028, 0.802805, 0.815276, 0.82743, 0.839255, 0.850737, 0.861867, 0.872632, 0.883022, 0.893027, 0.902635, 0.911838, 0.920627, 0.928992, 0.936925, 0.944418, 0.951463, 0.958054, 0.964184, 0.969846, 0.975036, 0.979746, 0.983974, 0.987715, 0.990964, 0.993719, 0.995977, 0.997736, 0.998993, 0.999748};
y = {-1.66533e-17, 3.65828e-05, 0.000146223, 0.000328595, 0.00058316, 0.00090917, 0.00130567, 0.00177151, 0.00230534, 0.00290565, 0.00357073, 0.00429874, 0.00508767, 0.00593537, 0.00683958, 0.00779793, 0.00880793, 0.00986703, 0.0109726, 0.0121219, 0.0133121, 0.0145406, 0.0158044, 0.0171006, 0.0184265, 0.0197789, 0.0211551, 0.0225521, 0.0239669, 0.0253966, 0.0268382, 0.0282887, 0.0297451, 0.0312044, 0.0326635, 0.0341194, 0.035569, 0.0370091, 0.0384366, 0.0398482, 0.0412407, 0.0426107, 0.0439549, 0.0452698, 0.0465521, 0.0477982, 0.0490045, 0.0501675, 0.0512835, 0.0523489, 0.0533601, 0.0543134, 0.0552053, 0.0560321, 0.0567903, 0.0574765, 0.0580873, 0.0586193, 0.0590696, 0.059435, 0.0597128, 0.0599003, 0.0599951, 0.0598982, 0.0597028, 0.0594075, 0.0590111, 0.0585128, 0.0579121, 0.0572087, 0.0564028, 0.0554948, 0.0544854, 0.0533756, 0.0521668, 0.0508605, 0.0494586, 0.0479633, 0.0463768, 0.0447018, 0.042941, 0.0410972, 0.0391735, 0.037173, 0.0350989, 0.0329544, 0.0307426, 0.0284669, 0.0261302, 0.0237357, 0.0212862, 0.0187844, 0.0162331, 0.0136345, 0.0109908, 0.008304, 0.0055757, 0.00280732, 0.0, -0.00280732, -0.0055757, -0.008304, -0.0109908, -0.0136345, -0.0162331, -0.0187844, -0.0212862, -0.0237357, -0.0261302, -0.0284669, -0.0307426, -0.0329544, -0.0350989, -0.037173, -0.0391735, -0.0410972, -0.042941, -0.0447018, -0.0463768, -0.0479633, -0.0494586, -0.0508605, -0.0521668, -0.0533756, -0.0544854, -0.0554948, -0.0564028, -0.0572087, -0.0579121, -0.0585128, -0.0590111, -0.0594075, -0.0597028, -0.0598982, -0.0599951, -0.0599003, -0.0597128, -0.059435, -0.0590696, -0.0586193, -0.0580873, -0.0574765, -0.0567903, -0.0560321, -0.0552053, -0.0543134, -0.0533601, -0.0523489, -0.0512835, -0.0501675, -0.0490045, -0.0477982, -0.0465521, -0.0452698, -0.0439549, -0.0426107, -0.0412407, -0.0398482, -0.0384366, -0.0370091, -0.035569, -0.0341194, -0.0326635, -0.0312044, -0.0297451, -0.0282887, -0.0268382, -0.0253966, -0.0239669, -0.0225521, -0.0211551, -0.0197789, -0.0184265, -0.0171006, -0.0158044, -0.0145406, -0.0133121, -0.0121219, -0.0109726, -0.00986703, -0.00880793, -0.00779793, -0.00683958, -0.00593537, -0.00508767, -0.00429874, -0.00357073, -0.00290565, -0.00230534, -0.00177151, -0.00130567, -0.00090917, -0.00058316, -0.000328595, -0.000146223, -3.65828e-05};

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
Line(400)  = {196, 402};                               // inlet front line
Line(401)  = {400, 2*point_id_top};                      // inlet top line
Line(402)  = {401, 2*point_id_bot+1};                      // inlet bottom line
Ellipse(403) = {400, 403, 403, 402};                   // inlet circle line top
Ellipse(404) = {402, 403, 403, 401};                   // inlet circle line bottom
Transfinite Curve {400, -401, -402} = n_around Using Progression progression_around;  // inlet lines points
Transfinite Curve {403} = n_inlet+1 Using Progression 1;        // inlet circle points
Transfinite Curve {404} = n_inlet Using Progression 1;        // inlet circle points
Curve Loop (500) = {-403, 401, 140:195, 400}; // inlet loop top
Curve Loop (501) = {-400, 196:252, -402, -404}; // inlet loop bottom
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
Mesh.SubdivisionAlgorithm = 1;
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
