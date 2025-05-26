import os.path

import numpy as np

nonUni = True
closed = True
ellipse= True

def rotate_points(content: str, aoa_deg: float):
    aoa_pi = aoa_deg/360*np.pi
    point_list = []
    x_coords = []
    y_coords = []
    lines_list = []
    for line in content.split("\n"):
        if "Point" in line:
            point_list.append(line.split('(')[1].split(')')[0])
            x_coords.append(float(line.split('{')[1].split(', ')[0]))
            y_coords.append(float(line.split('{')[1].split(', ')[1]))
        if "Line" in line:
            lines_list.append(line)
    if closed:
        x_coords = x_coords[:-1]
        y_coords = y_coords[:-1]
        point_list = point_list[:-1]
        lines_list = lines_list[:-2]
    x_coords = np.array(x_coords)
    y_coords = np.array(y_coords)
    x_new = x_coords*np.cos(aoa_pi) + y_coords*np.sin(aoa_pi)
    y_new = -x_coords*np.sin(aoa_pi) + y_coords*np.cos(aoa_pi)
    content_new = ""
    for _ in range(len(point_list)):
        content_new += f"Point({point_list[_]}) = " + '{' + f"{x_new[_]}, {y_new[_]}, 0, 1" + "};\n"
    for _ in lines_list:
        content_new += _ + "\n"
    if closed:
        content_new += "Line(198) = {198, 1};\n"
    return content_new


nth_point_to_top = 30  # this point will be connected to inlet circle
n_points_profile = 100
point_id_top     = n_points_profile - nth_point_to_top
point_id_bot     = n_points_profile + nth_point_to_top
n_channel        = int(1.5*n_points_profile)
inlet_center     = 0.3 if nonUni else 0.4

for deg in range(6):
    correction_bot = 1 if not closed else -1
    directoryname = f"{'nonUni' if nonUni else 'uni'}_{'closed' if closed else 'open'}/"
    filename = "NACA0012"
    # with open(filename + ".geo", 'r') as oldfile:
    content = open(directoryname + filename + ".geo", 'r').read()
    content = rotate_points(content, deg)
    # content = content.replace(", 1};", ", size_foil};")  # set element size to 1e-2 at airfoil
    setup_string = f"""
inlet_r      = 5;
inlet_front  = 3;
inlet_c      = {inlet_center};
outlet_c     = 6;
outlet_h     = inlet_r;
sponge_h     = inlet_r + 3;
size_foil    = 0.1;
size_in_out  = 1;
size_sponge  = 10;
point_id_top = {point_id_top};
point_id_bot = {point_id_bot};
n_around     = 40;
n_inlet      = 100 - point_id_top + 1;
n_sponge     = 9;
channel_l    = 1.5;
channel_h    = inlet_r;
n_channel    = {n_channel};
n_channel_top= n_channel-1;
n_channel_bot= n_channel+{correction_bot};
n_outlet_center = n_channel - {point_id_top};
progression_around = 1.3;
progression_sponge = 1.2;
    """
    f = open(directoryname + filename + f"_{deg}deg.geo", 'w')
    f.write(setup_string + '\n' + content + '\n')
    ### INLET
    f.write("\n/// INLET\n")
    f.write("Point(200) = {0, inlet_r, 0, size_in_out};       // inlet top\n")
    f.write("Point(201) = {0, -inlet_r, 0, size_in_out};      // inlet bottom\n")
    f.write("Point(202) = {inlet_c-inlet_front, 0, 0, size_in_out};     // inlet front\n")
    f.write("Point(203) = {inlet_c, 0, 0};                          // inlet center\n")
    f.write("Line(200)  = {100, 202};                               // inlet front line\n")
    f.write("Line(201)  = {200, point_id_top};                      // inlet top line\n")
    f.write("Line(202)  = {201, point_id_bot};                      // inlet bottom line\n")
    if ellipse:
        f.write("Ellipse(203) = {200, 203, 203, 202};                   // inlet circle line top\n")
        f.write("Ellipse(204) = {202, 203, 203, 201};                   // inlet circle line bottom\n")
    else:
        f.write("Circle(203) = {200, 203, 202};                   // inlet circle line top\n")
        f.write("Circle(204) = {202, 203, 201};                   // inlet circle line bottom\n")
    f.write("Transfinite Curve {200, -201, -202} = n_around Using Progression progression_around;  // inlet lines points\n")
    f.write("Transfinite Curve {203, 204} = n_inlet Using Progression 1;        // inlet circle points\n")
    transfinite_curve = "Transfinite Curve {"
    for point_id in np.arange(point_id_top, 100):
        transfinite_curve += f"{point_id}, "
    transfinite_curve = transfinite_curve.removesuffix(", ") + "} = 2 Using Progression 1; // inlet foil points top\n"
    f.write(transfinite_curve)
    transfinite_curve = "Transfinite Curve {"
    for point_id in np.arange(100, point_id_bot):
        transfinite_curve += f"{point_id}, "
    transfinite_curve = transfinite_curve.removesuffix(", ") + "} = 2 Using Progression 1; // inlet foil points bottom\n"
    f.write(transfinite_curve)
    curve_loop = "Curve Loop (200) = {-203, 201"
    for point_id in np.arange(point_id_top, n_points_profile):
        curve_loop += f", {point_id}"
    curve_loop += ", 200}; // inlet loop top\n"
    f.write(curve_loop)
    curve_loop = "Curve Loop (201) = {-200"
    for point_id in np.arange(n_points_profile, point_id_bot):
        curve_loop += f", {point_id}"
    curve_loop += ", -202, -204}; // inlet loop bottom\n"
    f.write(curve_loop)
    f.write("Plane Surface(1) = {200};                                // inlet surface top\n")
    f.write("Plane Surface(2) = {201};                                // inlet surface bottom\n")
    f.write("Transfinite Surface {1} = {202, 100, 70, 200};           // inlet surface top\n")
    f.write("Transfinite Surface {2} = {202, 201, 130, 100};          // inlet surface bottom\n")

    f.write("\n/// OUTLET\n")
    f.write("Point(210) = {outlet_c, outlet_h, 0, size_in_out};     // outlet top\n")
    f.write("Point(211) = {outlet_c, -outlet_h, 0, size_in_out};    // outlet bottom\n")
    f.write("Point(212) = {outlet_c, 0, 0, size_foil};              // outlet center\n")
    f.write("Line(210)  = {1, 212};                                 // outlet center line\n")
    f.write("Line(211)  = {210, 212};                               // outlet end top line\n")
    f.write("Line(212)  = {211, 212};                               // outlet end bottom line\n")
    # f.write("Line(213)  = {220, 210};                               // outlet top line\n")
    # f.write("Line(214)  = {221, 211};                               // outlet bottom line\n")
    f.write("Transfinite Curve {210, 213, 214} = n_outlet_center Using Progression 1;  // outlet center lines points\n")
    f.write("Transfinite Curve {-211, -212} = n_around Using Progression progression_around;              // outlet end lines points\n")

    f.write("\n/// CHANNEL\n")
    # f.write("Point(220) = {channel_l, channel_h, 0, size_in_out};   // channel end top\n")
    # f.write("Point(221) = {channel_l, -channel_h, 0, size_in_out};  // channel end bottom\n")
    f.write("Line(220)  = {200, 210};                               // channel top line\n")
    f.write("Line(221)  = {201, 211};                               // channel bottom line\n")
    # f.write("Line(222)  = {1, 220};                                 // channel end top line\n")
    # f.write("Line(223)  = {1, 221};                                 // channel end bottom line\n")
    # curve_loop = "Curve Loop(220) = {201"
    # for point_id in np.arange(-point_id_top+1, 0):
    #     curve_loop += f", {point_id}"
    # curve_loop += ", 222, -220}; // channel top\n"
    # f.write(curve_loop)
    # curve_loop = "Curve Loop(221) = {202"
    # for point_id in np.arange(point_id_bot, 200):
    #     curve_loop += f", {point_id}"
    # curve_loop += ", 223, -221}; // channel bottom\n"
    # f.write(curve_loop)
    f.write("Transfinite Curve {220} = n_channel_top Using Progression 1;      // channel top points\n")
    f.write("Transfinite Curve {221} = n_channel_bot Using Progression 1;      // channel bottom points\n")
    # f.write("Transfinite Curve {222, 223} = n_around Using Progression progression_around;      // channel end points\n")
    transfinite_curve = "Transfinite Curve {"
    for point_id in np.arange(1, point_id_top):
        transfinite_curve += f"{point_id}, "
    transfinite_curve = transfinite_curve.removesuffix(", ") + "} = 2 Using Progression 1; // channel foil points top\n"
    f.write(transfinite_curve)
    transfinite_curve = "Transfinite Curve {"
    for point_id in np.arange(point_id_bot, 2*n_points_profile-[1 if closed else 0][0]):
        transfinite_curve += f"{point_id}, "
    transfinite_curve = transfinite_curve.removesuffix(", ") + "} = 2 Using Progression 1; // channel foil points bottom\n"
    f.write(transfinite_curve)
    # f.write("Plane Surface(3) = {220};                              // channel top\n")
    # f.write("Plane Surface(4) = {221};                              // channel bottom\n")
    # f.write("Transfinite Surface {3} = {1, 220, 200, point_id_top}; // channel surface top\n")
    # f.write("Transfinite Surface {4} = {point_id_bot, 201, 221, 1}; // channel surface bottom\n")

    f.write("\n/// OUTLET / CHANNEL SURFACES\n")
    curve_loop = "Curve Loop(222) = {-201, 220, 211, -210"
    for line_id in np.arange(1, point_id_top):
        curve_loop += f", {line_id}"
    curve_loop += "}; // channel top\n"
    f.write(curve_loop)
    curve_loop = "Curve Loop(223) = {202"
    for line_id in np.arange(point_id_bot, 200-[1 if closed else 0][0]):
        curve_loop += f", {line_id}"
    curve_loop += ", 210, -212, -221}; // channel bottom\n"
    f.write(curve_loop)
    f.write("Plane Surface(5) = {222};                               // channel top\n")
    f.write("Plane Surface(6) = {223};                               // channel bottom\n")
    f.write("Transfinite Surface {5} = {point_id_top, 212, 210, 200}; // channel surface top\n")
    f.write("Transfinite Surface {6} = {point_id_bot, 201, 211, 212}; // channel surface bottom\n")

    f.write("\n/// SPONGE\n")
    f.write("Point(230) = {outlet_c, sponge_h, 0, size_sponge};      // sponge top back\n")
    f.write("Point(231) = {outlet_c, -sponge_h, 0, size_sponge};     // sponge bottom back\n")
    f.write("Point(232) = {outlet_c/2, sponge_h, 0, size_sponge};    // sponge top front\n")
    f.write("Point(233) = {outlet_c/2, -sponge_h, 0, size_sponge};   // sponge bottom front\n")
    f.write("Line(230)  = {210, 230};                                // sponge outlet top line\n")
    f.write("Line(231)  = {211, 231};                                // sponge outlet bottom line\n")
    f.write("Line(232)  = {200, 232};                                // sponge top line\n")
    f.write("Line(233)  = {201, 233};                                // sponge bottom line\n")
    f.write("Line(234)  = {232, 230};                                // sponge front top line\n")
    f.write("Line(235)  = {233, 231};                                // sponge front bottom line\n")
    f.write("Curve Loop(230) = {232, 234, -230, -220};               // sponge surface top\n")
    f.write("Curve Loop(231) = {221, 231, -235, -233};               // sponge surface top\n")
    f.write("Transfinite Curve {230, 231, 232, 233} = n_sponge*2 Using Progression progression_sponge;  // sponge lines points\n")
    f.write("Transfinite Curve {234, 235} = n_sponge Using Progression 1; // sponge top lines points\n")
    f.write("Plane Surface(7) = {230};                               // sponge top\n")
    f.write("Plane Surface(8) = {231};                               // sponge bottom\n")
    # f.write("Transfinite Surface {7} Right;                           // sponge surface top\n")
    # f.write("Transfinite Surface {8} Right;                           // sponge surface bottom\n")

    f.write("\n/// MESH SIZES\n")
    # meshsize_airfoil = "MeshSize {"
    # for point_id in range(200):
    #     meshsize_airfoil += f"{point_id}, "
    # meshsize_airfoil = meshsize_airfoil.removesuffix(", ") + "} = size_foil; // mesh size airfoil nodes\n"
    # f.write(meshsize_airfoil)
    f.write("Mesh.Algorithm = 6;\n")
    f.write("Mesh.SubdivisionAlgorithm = 0;\n")
    f.write("Mesh.RecombineAll = 1;\n")
    # f.write("Mesh.QualityFactor = 0.00001;\n")
    # f.write("Mesh.Optimize = 1;\n")
    # f.write("Mesh.OptimizeNetgen = 1;\n")
    # meshsize_around = "MeshSize {"
    # for point_id in np.arange(200, 240):
    #     meshsize_around += f"{point_id}, "
    # meshsize_around = meshsize_around.removesuffix(", ") + "} = size_around; // mesh size airfoil nodes\n"
    # f.write(meshsize_around)

    f.write("\n/// BOUNDARIES\n")
    f.write('Physical Curve("Inlet", 300) = {233, 204, 203, 232};\n')
    f.write('Physical Curve("Sponge", 301) = {234, 235};\n')
    f.write('Physical Curve("Outlet", 302) = {231, 212, 211, 230};\n')
    f.write('Physical Curve("BB_BC", 303) = {1:' + str(198 if closed else 199) + '};\n')
    f.write("Physical Surface(236) = {1, 2, 6, 5, 7, 8};\n")

    # f.write('Mesh.MeshSizeMin = 0\n')
    # f.write('Mesh.MeshSizeMax = 1e+22;\n')
    # f.write("Mesh 1;\n")
    f.write("Mesh 2;\n")
    f.write(f'Save "{filename+f"_{deg}deg"}.msh";')
    # f.write("Exit;\n")
