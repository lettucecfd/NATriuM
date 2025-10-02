//FILE newNonUni 10deg res40

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
point_id_front = 40;
point_id_bot = 53;
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

x={0.9961946980917455, 0.9946790849483617, 0.9901397792970809, 0.982605708265676, 0.9721209536954039, 0.9587500817374683, 0.9425719042176021, 0.9236858915829695, 0.9022052866950232, 0.8782618243357653, 0.8519991064678327, 0.823578840869913, 0.793173344839402, 0.7609668900100729, 0.7271578690529573, 0.6919523827294435, 0.6555634180080634, 0.6182158290379837, 0.5801341213453998, 0.5415529367149734, 0.5027045759189995, 0.4638262335475204, 0.42515527850206397, 0.3869242730221532, 0.3493680332039356, 0.3127117618195485, 0.2771816203542996, 0.24299233889157595, 0.21035263286504824, 0.17946429401971534, 0.15051621324371578, 0.12368619864623859, 0.09914205890796388, 0.0770344556066909, 0.05750279322978377, 0.040667548565342874, 0.02663560029866631, 0.015492819870649469, 0.007309438555774413, 0.002134506072635025, 0.0, 0.0009357659868837347, 0.0049557105671311575, 0.012032039637625459, 0.022122153004736085, 0.035162791853400786, 0.05107645169402794, 0.06976479510410871, 0.09111518931239004, 0.11499607384835508, 0.1412632366591679, 0.16975374978574223, 0.20029276841414256, 0.2326926218366287, 0.2667506467507744, 0.30225516442765543, 0.3389865640622911, 0.37671275757886663, 0.4152006981883975, 0.45420702853194683, 0.49349012217274607, 0.5328031973890506, 0.57189929814763, 0.6105365365444876, 0.6484713809092009, 0.6854700872068437, 0.7212992600254596, 0.7557375454452134, 0.7885706500648981, 0.8195925115081206, 0.8486108398127747, 0.875445299353132, 0.8999268612681139, 0.921904253889722, 0.9412382470420774, 0.957808974027279, 0.9715106891846848, 0.9822588284095403, 0.9899844677635045, 0.9946400391756106};
y={-0.08682408883346517, -0.08646884683186289, -0.08540875841263966, -0.08365730288278113, -0.08123824612186553, -0.07818212030820969, -0.07452866877005968, -0.07032242393725083, -0.06561111097401143, -0.06044904769513049, -0.0548925824030377, -0.0489976490366125, -0.04282509099400995, -0.03643679356245103, -0.02989388036925128, -0.023261135430160652, -0.01660493912812551, -0.009993702333732882, -0.003499790311868131, 0.002805540628844236, 0.008847253433658935, 0.01454888832881198, 0.019835962737688555, 0.0246364054012672, 0.028880937208072575, 0.03250510168181611, 0.035455316930529765, 0.0376849745427591, 0.039158948459992085, 0.03985467039932256, 0.039762650797983176, 0.038884327968018564, 0.03723296987077111, 0.034832304688353605, 0.031715011280541563, 0.0279154181647984, 0.023473023790129558, 0.018428161312831384, 0.01281455746193595, 0.006664796476617042, 0.0, -0.006932388318401782, -0.013883535643653574, -0.020827110887300024, -0.027722541993994678, -0.03452446780680177, -0.04117831601885627, -0.04762670241885303, -0.05381498165254163, -0.05968685870798382, -0.06519290547277212, -0.07029102838807144, -0.0749490534143901, -0.07914365060891637, -0.08286213767449008, -0.08610300090413549, -0.08887499340657568, -0.09119175799291526, -0.09307781205857943, -0.09456075915240347, -0.09567134226712411, -0.09644184747218308, -0.09690653803417133, -0.0970991227415494, -0.09704918234030172, -0.09678914301445031, -0.09634747655371872, -0.095752708038322, -0.09503298171852237, -0.09421417064156895, -0.09332534058910369, -0.09239659923183459, -0.0914550549111484, -0.09053135599918008, -0.08965620415855596, -0.08885700771671727, -0.08816041334119967, -0.08759192520968057, -0.08717044107257303, -0.08691173899328271};
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
Save "NACA0012_10deg.msh";