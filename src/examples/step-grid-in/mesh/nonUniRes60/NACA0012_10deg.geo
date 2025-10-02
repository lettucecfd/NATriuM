//FILE newNonUni 10deg res60

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

x={0.9961946980917455, 0.9935007216435431, 0.9901397792970809, 0.9854477194267082, 0.9794359734343455, 0.9721209536954039, 0.9635228830526015, 0.9536648857770073, 0.9425719042176021, 0.930275941241995, 0.9168077420558693, 0.9022052866950232, 0.8865081989616799, 0.8697572731642464, 0.8519991064678327, 0.8332815536994748, 0.8136550224915897, 0.793173344839402, 0.771891000309864, 0.7498656313655071, 0.7271578690529573, 0.703829253457796, 0.6799433170550008, 0.6555634180080634, 0.6307581569211947, 0.6055919753083276, 0.5801341213453998, 0.5544541046755772, 0.5286210863190551, 0.5027045759189995, 0.47677464353456067, 0.45090179500360295, 0.42515527850206397, 0.3996045165173671, 0.37431893153693563, 0.3493680332039356, 0.3248164373437854, 0.3007330060279239, 0.2771816203542996, 0.25422633573234654, 0.2319295638050734, 0.21035263286504824, 0.18955431839945844, 0.1695904073113971, 0.15051621324371578, 0.1323837997879435, 0.115242939197457, 0.09914205890796388, 0.08412404496569696, 0.0702309615336765, 0.05750279322978377, 0.04597255130901344, 0.03567441753254376, 0.02663560029866631, 0.018881228460533774, 0.012433180819976297, 0.007309438555774413, 0.003522914717475448, 0.0010842684980022974, 0.0, 0.0002805182383833939, 0.0019342398386711346, 0.0049557105671311575, 0.009335665722724526, 0.015063109682244367, 0.022122153004736085, 0.030492834314709987, 0.04015246511981033, 0.05107645169402794, 0.06323322694846703, 0.07659204567744433, 0.09111518931239004, 0.10676303405108085, 0.12349478558510985, 0.1412632366591679, 0.16002053244621875, 0.1797151322891898, 0.20029276841414256, 0.22169765391798385, 0.24387101331352626, 0.2667506467507744, 0.2902734453407774, 0.31437395918221744, 0.3389865640622911, 0.36404065479470427, 0.3894693495047274, 0.4152006981883975, 0.44116266368600865, 0.46728320883785524, 0.49349012217274607, 0.51971045749202, 0.5458708328183023, 0.57189929814763, 0.597723554853069, 0.6232716529306565, 0.6484713809092009, 0.6732556826024875, 0.6975536913569936, 0.7212992600254596, 0.7444264157721113, 0.7668711781505699, 0.7885706500648981, 0.8094649230032533, 0.8294969027264004, 0.8486108398127747, 0.8667535376461912, 0.8838752239732732, 0.8999268612681139, 0.9148655634844804, 0.9286492664593526, 0.9412382470420774, 0.95259949397766, 0.962699261283636, 0.9715106891846848, 0.9790090846063676, 0.9851728302140821, 0.9899844677635045, 0.9934315199838014};
y={-0.08682408883346517, -0.08619380402003496, -0.08540875841263966, -0.08431642888672991, -0.08292378892281929, -0.08123824612186553, -0.07926953257833097, -0.07702862945751891, -0.07452866877005968, -0.07178240169663071, -0.06880545838954447, -0.06561111097401143, -0.06221675961771812, -0.05863791411529075, -0.0548925824030377, -0.050995893445632344, -0.04696817967288763, -0.04282509099400995, -0.03858631853659295, -0.034269836703347586, -0.02989388036925128, -0.025477759831944946, -0.021040959037247113, -0.01660493912812551, -0.012188629572268274, -0.00781258973815737, -0.003499790311868131, 0.0007297638162479039, 0.004852113361540768, 0.008847253433658935, 0.012690149325289043, 0.016360709321990692, 0.019835962737688555, 0.023094916083456485, 0.026116575870368343, 0.028880937208072575, 0.03136941792523641, 0.03356605434989553, 0.035455316930529765, 0.03702365331276816, 0.038261639184715826, 0.039158948459992085, 0.03971045851735579, 0.039911307207832224, 0.039762650797983176, 0.039262940230913077, 0.038418795710589655, 0.03723296987077111, 0.03571367928262161, 0.033871204538543344, 0.031715011280541563, 0.029255172919243403, 0.026503653238336297, 0.023473023790129558, 0.02017504117653555, 0.016618756676009348, 0.01281455746193595, 0.008772276229569942, 0.004498953526619761, 0.0, -0.004617902528321608, -0.009247898588199664, -0.013883535643653574, -0.01851603666519823, -0.023133485179447043, -0.027722541993994678, -0.032270509218655055, -0.0367614626952518, -0.04117831601885627, -0.045503375015918, -0.04972100953412455, -0.05381498165254163, -0.057767891203482705, -0.0615642283923216, -0.06519290547277212, -0.06863879348016085, -0.07189441176614465, -0.0749490534143901, -0.07779785627592253, -0.08043569772950074, -0.08286213767449008, -0.0850756605875924, -0.08707853169163034, -0.08887499340657568, -0.09046885403195605, -0.09186723495937824, -0.09307781205857943, -0.09410924979787141, -0.09497021264556597, -0.09567134226712411, -0.09622022770819416, -0.09662940100729753, -0.09690653803417133, -0.09706326905285123, -0.09710726993307436, -0.09704918234030172, -0.0968981048632893, -0.09666081159728851, -0.09634747655371872, -0.09596629654685018, -0.09552564203913069, -0.09503298171852237, -0.09449604474525385, -0.09392479794896935, -0.09332534058910369, -0.09270897539023105, -0.09208232255631917, -0.0914550549111484, -0.09083538902588294, -0.09023360549292497, -0.08965620415855596, -0.08911325843340302, -0.08861178910828063, -0.08816041334119967, -0.08776594474119957, -0.08743446879101154, -0.08717044107257303, -0.08697875128826565};
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