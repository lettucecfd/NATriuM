import numpy as np
from glob import glob as glob
import os

##### CAREFUL! For low resolutions, this is not valid! Had to change front and bot point_ids for res <= 50  

soureceDir = "/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/mesh/varyRefinement"
aoa_deg = 0

soureceFile = soureceDir + "/naca0012_res100_combined.txt"

res = 100
targetDir = f"/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/mesh/nonUniRectangular"
if not os.path.exists(targetDir):
  os.mkdir(targetDir)

x = []
y = []
with open(soureceFile, "r") as file:
  string = file.read()
  string = string.splitlines()[1:-1]
  for line in string:
    if line not in ['', ' ']:
      line = line.split(' ')
      line = [l for l in line if l != '']
      x.append(float(line[0]))
      y.append(float(line[1]))
tol = 1e-5/res  # TODO may need to adapt depending on resolution
for i in range(len(x)):
  if abs(float(y[i]) - max(y)) < tol:
    point_id_top = i+2
  if abs(float(x[i])) < tol:
    point_id_front = i+1  # TODO: may need to adapt depending on resolution
n_channel = int(res/100*400)
point_id_bot = int(2*point_id_front - point_id_top)
aoa_pi = aoa_deg/360*np.pi
x = np.array(x)
y = np.array(y)
x = x*np.cos(aoa_pi) + y*np.sin(aoa_pi)
y = -x*np.sin(aoa_pi) + y*np.cos(aoa_pi)

xstring = "x={" + ''.join([str(xi) + ", " for xi in x[:-1]]) + str(x[-1]) + "};"
ystring = "y={" + ''.join([str(yi) + ", " for yi in y[:-1]]) + str(y[-1]) + "};"

header = f"//FILE newNonUni {aoa_deg}deg res{res}"
variables = f"""

inlet_r      = 4;
inlet_front  = 3;
inlet_c      = 0.3;
outlet_c     = 6;
outlet_h     = inlet_r;
sponge_h     = inlet_r * 3;
size_foil    = 0.1;
size_in_out  = 1;
size_sponge  = 1;
point_id_top = {point_id_top};
point_id_front = {point_id_front};
point_id_bot = {point_id_bot};
n_around     = {res/100*70};
n_foil       = 2;  // number of points per foil section; 2 for just the section, 3 to split once, ...
n_inlet      = point_id_front - point_id_top + 1;  // top and bottom, each; 31 default points; additional 31-1 (30 segments) points for each n_foil above 2
channel_l    = 1.5;
channel_h    = inlet_r;
n_channel    = {n_channel};  // number of points on top and bottom
n_outlet_center = n_channel - point_id_top + 1;
progression_around = 1;

"""
domain = """
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
"""
saveline = f'Save "NACA0012_{aoa_deg}deg.msh";'

with open(targetDir + f"/NACA0012_{aoa_deg}deg.geo", 'w') as file:
  file.write(header)
  file.write(variables)
  file.write(xstring)
  file.write("\n")
  file.write(ystring)
  file.write(domain)
  file.write(saveline)
