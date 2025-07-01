//FILE newNonUni 10deg res50

inlet_r      = 4;
inlet_front  = 3;
inlet_c      = 0.3;
outlet_c     = 6;
outlet_h     = inlet_r;
sponge_h     = inlet_r * 3;
size_foil    = 0.1;
size_in_out  = 1;
size_sponge  = 1;
point_id_top = 34;
point_id_front = 50;
point_id_bot = 66;
n_around     = 50;
n_foil       = 2;  // number of points per foil section; 2 for just the section, 3 to split once, ...
n_inlet      = point_id_front - point_id_top + 1;  // top and bottom, each; 31 default points; additional 31-1 (30 segments) points for each n_foil above 2
channel_l    = 1.5;
channel_h    = inlet_r;
n_channel    = 175;  // number of points on top and bottom
n_outlet_center = n_channel - point_id_top + 1;
progression_around = 1.05;
progression_sponge_front = 1.05;
progression_sponge_back = 1.05;
n_sponge_front= 57;
n_sponge_back= 50;

x={0.9961946980917455, 0.9952239171959419, 0.9923165554820214, 0.9874839571961621, 0.98074344375273, 0.9721209536954039, 0.961651217008659, 0.9493745922221886, 0.9353363316834333, 0.9195935549204687, 0.9022052866950232, 0.8832412484432762, 0.8627724192633337, 0.8408812968109645, 0.8176507273649086, 0.793173344839402, 0.7675426921876352, 0.7408613652708248, 0.7132318727944444, 0.6847617915932008, 0.6555634180080634, 0.6257512150808828, 0.5954402934252367, 0.564749013706026, 0.5337963841598777, 0.5027045759189995, 0.47159350625897395, 0.44058352823409785, 0.40979655150935, 0.37935193294432923, 0.3493680332039356, 0.31996048583620373, 0.2912442768174408, 0.2633326612017862, 0.23633408738137435, 0.21035263286504824, 0.18549272378462212, 0.16185016914259823, 0.13952039183698003, 0.11858902713027558, 0.09914205890796388, 0.08125395249786056, 0.06499843850031144, 0.05043820221824055, 0.03763428518814676, 0.02663560029866631, 0.017486166620687756, 0.010223590278197345, 0.0048765987986412644, 0.0014644310227263735, 0.0, 0.0005020573113067321, 0.0029793925905102414, 0.00742101021440365, 0.013810286014562525, 0.022122153004736085, 0.03232249129064798, 0.04437163958854706, 0.05822486848544894, 0.07382566854347765, 0.09111518931239004, 0.11002370935600672, 0.13048024957622192, 0.1524014407598391, 0.17570356507069063, 0.20029276841414256, 0.22607359041066355, 0.2529412563054685, 0.28079064272080123, 0.30951103521521606, 0.3389865640622911, 0.3690991056789802, 0.399729191664568, 0.43075410356701693, 0.46205012673959084, 0.49349012217274607, 0.5249493790250485, 0.5563027506763505, 0.5874233595843361, 0.6181871424792991, 0.6484713809092009, 0.6781560835388704, 0.7071226038508045, 0.7352541134054115, 0.7624390262038179, 0.7885706500648981, 0.8135423799632696, 0.8372564894700893, 0.8596163355069555, 0.8805354112539324, 0.8999268612681139, 0.9177162202216841, 0.9338297575142973, 0.948204962154515, 0.9607814026960374, 0.9715106891846848, 0.9803494997955107, 0.987260838494728, 0.9922168493123181, 0.995198990653516};
y={-0.08682408883346517, -0.08659702386161905, -0.08591626306652485, -0.08478976852727076, -0.08322602326747795, -0.08123824612186553, -0.07884241453910126, -0.07605898130668433, -0.07290982193113263, -0.06941884240660473, -0.06561111097401143, -0.06151561673496445, -0.05715846981944064, -0.05257176523935576, -0.04778364361232722, -0.04282509099400995, -0.03772816886272229, -0.032524731515007545, -0.027247621845983554, -0.021929128272637738, -0.01660493912812551, -0.01130893919663073, -0.006079054480724468, -0.0009505056595205535, 0.0040374453694800375, 0.008847253433658935, 0.01343858121285172, 0.01777603437976691, 0.021819101966062404, 0.025532476468535556, 0.028880937208072575, 0.03183248977354997, 0.03435918097223145, 0.03643619705528221, 0.03804113559146079, 0.039158948459992085, 0.0397785419343994, 0.03989355809930391, 0.0395046460534378, 0.03861530033357084, 0.03723296987077111, 0.035371088478045765, 0.03304614258821452, 0.03027458510049968, 0.027076996956599642, 0.023473023790129558, 0.01948373395372338, 0.015125773751675058, 0.010417559503652642, 0.005372357354676097, 0.0, -0.005543748106033357, -0.011102254268193348, -0.016663602013093394, -0.022211399528515523, -0.027722541993994678, -0.0331741317708409, -0.038537807282958234, -0.04378576143238116, -0.04888716803492896, -0.05381498165254163, -0.05854021318360743, -0.06303675249814954, -0.06728239056999816, -0.07125887276728186, -0.0749490534143901, -0.07834262020164565, -0.0814324337968702, -0.08421530926218382, -0.08669454772928462, -0.08887499340657568, -0.09076428265000627, -0.09237391938573519, -0.09371811286337856, -0.0948109855084636, -0.09567134226712411, -0.09631321874079848, -0.09675559352379814, -0.09701430576653305, -0.0971074322888289, -0.09704918234030172, -0.09685699143855796, -0.0965444275309944, -0.09612720941033479, -0.09561852419402317, -0.09503298171852237, -0.09438420322172064, -0.09368757995688032, -0.09295760140277798, -0.0922076480819293, -0.0914550549111484, -0.09071325570344241, -0.08999873689163106, -0.08932597417778747, -0.08870862831358782, -0.08816041334119967, -0.08769448882466023, -0.08732058087824125, -0.08704721983586478, -0.08687976305395403};
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