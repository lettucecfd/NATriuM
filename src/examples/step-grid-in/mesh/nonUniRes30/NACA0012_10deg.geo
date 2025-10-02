//FILE newNonUni 10deg res30

inlet_r      = 4;
inlet_front  = 3;
inlet_c      = 0.3;
outlet_c     = 6;
outlet_h     = inlet_r;
sponge_h     = inlet_r * 3;
size_foil    = 0.1;
size_in_out  = 1;
size_sponge  = 1;
point_id_top = 21;
point_id_front = 30;
point_id_bot = 39;
n_around     = 30;
n_foil       = 2;  // number of points per foil section; 2 for just the section, 3 to split once, ...
n_inlet      = point_id_front - point_id_top + 1;  // top and bottom, each; 31 default points; additional 31-1 (30 segments) points for each n_foil above 2
channel_l    = 1.5;
channel_h    = inlet_r;
n_channel    = 105;  // number of points on top and bottom
n_outlet_center = n_channel - point_id_top + 1;
progression_around = 1.05;
progression_sponge_front = 1.05;
progression_sponge_back = 1.05;
n_sponge_front= 57;
n_sponge_back= 50;

x={0.9961946980917455, 0.9935007216435431, 0.9854477194267082, 0.9721209536954039, 0.9536648857770073, 0.930275941241995, 0.9022052866950232, 0.8697572731642464, 0.8332815536994748, 0.793173344839402, 0.7498656313655071, 0.703829253457796, 0.6555634180080634, 0.6055919753083276, 0.5544541046755772, 0.5027045759189995, 0.45090179500360295, 0.3996045165173671, 0.3493680332039356, 0.3007330060279239, 0.25422633573234654, 0.21035263286504824, 0.1695904073113971, 0.1323837997879435, 0.09914205890796388, 0.0702309615336765, 0.04597255130901344, 0.02663560029866631, 0.012433180819976297, 0.003522914717475448, 0.0, 0.0019342398386711346, 0.009335665722724526, 0.022122153004736085, 0.04015246511981033, 0.06323322694846703, 0.09111518931239004, 0.12349478558510985, 0.16002053244621875, 0.20029276841414256, 0.24387101331352626, 0.2902734453407774, 0.3389865640622911, 0.3894693495047274, 0.44116266368600865, 0.49349012217274607, 0.5458708328183023, 0.597723554853069, 0.6484713809092009, 0.6975536913569936, 0.7444264157721113, 0.7885706500648981, 0.8294969027264004, 0.8667535376461912, 0.8999268612681139, 0.9286492664593526, 0.95259949397766, 0.9715106891846848, 0.9851728302140821, 0.9934315199838014};
y={-0.08682408883346517, -0.08619380402003496, -0.08431642888672991, -0.08123824612186553, -0.07702862945751891, -0.07178240169663071, -0.06561111097401143, -0.05863791411529075, -0.050995893445632344, -0.04282509099400995, -0.034269836703347586, -0.025477759831944946, -0.01660493912812551, -0.00781258973815737, 0.0007297638162479039, 0.008847253433658935, 0.016360709321990692, 0.023094916083456485, 0.028880937208072575, 0.03356605434989553, 0.03702365331276816, 0.039158948459992085, 0.039911307207832224, 0.039262940230913077, 0.03723296987077111, 0.033871204538543344, 0.029255172919243403, 0.023473023790129558, 0.016618756676009348, 0.008772276229569942, 0.0, -0.009247898588199664, -0.01851603666519823, -0.027722541993994678, -0.0367614626952518, -0.045503375015918, -0.05381498165254163, -0.0615642283923216, -0.06863879348016085, -0.0749490534143901, -0.08043569772950074, -0.0850756605875924, -0.08887499340657568, -0.09186723495937824, -0.09410924979787141, -0.09567134226712411, -0.09662940100729753, -0.09706326905285123, -0.09704918234030172, -0.09666081159728851, -0.09596629654685018, -0.09503298171852237, -0.09392479794896935, -0.09270897539023105, -0.0914550549111484, -0.09023360549292497, -0.08911325843340302, -0.08816041334119967, -0.08743446879101154, -0.08697875128826565};
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