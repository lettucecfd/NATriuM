import gmsh
import numpy as np
import meshio
import os
import shapely

gmsh.initialize()
gmsh.model.add("clipped_cartesian")

# Parameters
Lx = 9
Ly = 7
xmin = -1.5
ymin = -Ly/2
xmax = xmin + Lx
ymax = ymin + Ly
dx_foil = 0.001
dx_around = 0.35
prog_y = 1.3
prog_x = 1.2

# foil coordinates
x = [1.0, 0.999748, 0.998993, 0.997736, 0.995977, 0.993719, 0.990964, 0.987715, 0.983974, 0.979746, 0.975036, 0.969846, 0.964184, 0.958054, 0.951463, 0.944418, 0.936925, 0.928992, 0.920627, 0.911838, 0.902635, 0.893027, 0.883022, 0.872632, 0.861867, 0.850737, 0.839255, 0.82743, 0.815276, 0.802805, 0.790028, 0.77696, 0.763613, 0.75, 0.736136, 0.722033, 0.707708, 0.693173, 0.678443, 0.663534, 0.64846, 0.633237, 0.617879, 0.602403, 0.586824, 0.571157, 0.555419, 0.539625, 0.523791, 0.507933, 0.492067, 0.476209, 0.460375, 0.444581, 0.428843, 0.413176, 0.397597, 0.382121, 0.366763, 0.35154, 0.336466, 0.321557, 0.292292, 0.277967, 0.263864, 0.25, 0.236387, 0.22304, 0.209972, 0.197195, 0.184724, 0.17257, 0.160745, 0.149263, 0.138133, 0.127368, 0.116978, 0.106973, 0.0973649, 0.0881617, 0.0793732, 0.0710083, 0.0630753, 0.0555823, 0.0485367, 0.0419458, 0.035816, 0.0301537, 0.0249644, 0.0202535, 0.0160256, 0.0122851, 0.00903565, 0.00628056, 0.00402259, 0.00226404, 0.00100666, 0.000251729, 0.0, 0.000251729, 0.00100666, 0.00226404, 0.00402259, 0.00628056, 0.00903565, 0.0122851, 0.0160256, 0.0202535, 0.0249644, 0.0301537, 0.035816, 0.0419458, 0.0485367, 0.0555823, 0.0630753, 0.0710083, 0.0793732, 0.0881617, 0.0973649, 0.106973, 0.116978, 0.127368, 0.138133, 0.149263, 0.160745, 0.17257, 0.184724, 0.197195, 0.209972, 0.22304, 0.236387, 0.25, 0.263864, 0.277967, 0.306827, 0.321557, 0.336466, 0.35154, 0.366763, 0.382121, 0.397597, 0.413176, 0.428843, 0.444581, 0.460375, 0.476209, 0.492067, 0.507933, 0.523791, 0.539625, 0.555419, 0.571157, 0.586824, 0.602403, 0.617879, 0.633237, 0.64846, 0.663534, 0.678443, 0.693173, 0.707708, 0.722033, 0.736136, 0.75, 0.763613, 0.77696, 0.790028, 0.802805, 0.815276, 0.82743, 0.839255, 0.850737, 0.861867, 0.872632, 0.883022, 0.893027, 0.902635, 0.911838, 0.920627, 0.928992, 0.936925, 0.944418, 0.951463, 0.958054, 0.964184, 0.969846, 0.975036, 0.979746, 0.983974, 0.987715, 0.990964, 0.993719, 0.995977, 0.997736, 0.998993, 0.999748];
y = [-1.66533e-17, 3.65828e-05, 0.000146223, 0.000328595, 0.00058316, 0.00090917, 0.00130567, 0.00177151, 0.00230534, 0.00290565, 0.00357073, 0.00429874, 0.00508767, 0.00593537, 0.00683958, 0.00779793, 0.00880793, 0.00986703, 0.0109726, 0.0121219, 0.0133121, 0.0145406, 0.0158044, 0.0171006, 0.0184265, 0.0197789, 0.0211551, 0.0225521, 0.0239669, 0.0253966, 0.0268382, 0.0282887, 0.0297451, 0.0312044, 0.0326635, 0.0341194, 0.035569, 0.0370091, 0.0384366, 0.0398482, 0.0412407, 0.0426107, 0.0439549, 0.0452698, 0.0465521, 0.0477982, 0.0490045, 0.0501675, 0.0512835, 0.0523489, 0.0533601, 0.0543134, 0.0552053, 0.0560321, 0.0567903, 0.0574765, 0.0580873, 0.0586193, 0.0590696, 0.059435, 0.0597128, 0.0599003, 0.0599951, 0.0598982, 0.0597028, 0.0594075, 0.0590111, 0.0585128, 0.0579121, 0.0572087, 0.0564028, 0.0554948, 0.0544854, 0.0533756, 0.0521668, 0.0508605, 0.0494586, 0.0479633, 0.0463768, 0.0447018, 0.042941, 0.0410972, 0.0391735, 0.037173, 0.0350989, 0.0329544, 0.0307426, 0.0284669, 0.0261302, 0.0237357, 0.0212862, 0.0187844, 0.0162331, 0.0136345, 0.0109908, 0.008304, 0.0055757, 0.00280732, 0.0, -0.00280732, -0.0055757, -0.008304, -0.0109908, -0.0136345, -0.0162331, -0.0187844, -0.0212862, -0.0237357, -0.0261302, -0.0284669, -0.0307426, -0.0329544, -0.0350989, -0.037173, -0.0391735, -0.0410972, -0.042941, -0.0447018, -0.0463768, -0.0479633, -0.0494586, -0.0508605, -0.0521668, -0.0533756, -0.0544854, -0.0554948, -0.0564028, -0.0572087, -0.0579121, -0.0585128, -0.0590111, -0.0594075, -0.0597028, -0.0598982, -0.0599951, -0.0599003, -0.0597128, -0.059435, -0.0590696, -0.0586193, -0.0580873, -0.0574765, -0.0567903, -0.0560321, -0.0552053, -0.0543134, -0.0533601, -0.0523489, -0.0512835, -0.0501675, -0.0490045, -0.0477982, -0.0465521, -0.0452698, -0.0439549, -0.0426107, -0.0412407, -0.0398482, -0.0384366, -0.0370091, -0.035569, -0.0341194, -0.0326635, -0.0312044, -0.0297451, -0.0282887, -0.0268382, -0.0253966, -0.0239669, -0.0225521, -0.0211551, -0.0197789, -0.0184265, -0.0171006, -0.0158044, -0.0145406, -0.0133121, -0.0121219, -0.0109726, -0.00986703, -0.00880793, -0.00779793, -0.00683958, -0.00593537, -0.00508767, -0.00429874, -0.00357073, -0.00290565, -0.00230534, -0.00177151, -0.00130567, -0.00090917, -0.00058316, -0.000328595, -0.000146223, -3.65828e-05];

foil_top = max(y)+5*dx_foil
foil_bottom = min(y)-5*dx_foil
foil_front = min(x)-5*dx_foil
foil_back = max(x)+5*dx_foil

# Transfinite (structured) meshing
nx_foil = int((max(x)-min(x)) / dx_foil) + 1
ny_foil = int((max(y)-min(y)) / dx_foil) + 1
nx_around = int(Lx / dx_around) + 1
ny_around = int(Lx / dx_around) + 1


# Geometry points (as in your original setup)
p = {}  # map labels to point tags for reuse

def addp(label, x, y, dx, i):
    p[label] = gmsh.model.geo.addPoint(x, y, 0, meshSize=dx, tag=i)

# Corner + intermediate points

xcorners = [xmin, foil_front, foil_back, xmax]
ycorners = [ymin, foil_bottom, foil_top, ymax]
xlabels = [0, 1, 2, 3]
ylabels = [0, 1, 2, 3]
size = [dx_around, dx_foil, dx_foil, dx_around]
i = 0
for xi, ix in zip(xcorners, xlabels):
    for yi, iy in zip(ycorners, ylabels):
        addp(f"x{ix}_y{iy}", xi, yi, size[ix]*size[iy], i)
        i += 1

l = {}

# Lines (in counter-clockwise loops)
def add_line(label, start, end, n, prog=False):
    l[label] = gmsh.model.geo.addLine(p[start], p[end])

# bottom left
add_line("l1", "x0_y0", "x0_y1", ny_around, True)
add_line("l2", "x0_y1", "x1_y1", nx_around, True)
add_line("l3", "x1_y1", "x1_y0", ny_around, True)
add_line("l4", "x1_y0", "x0_y0", nx_around, True)

# mid left
add_line("l5", "x0_y1", "x0_y2", ny_foil)
add_line("l6", "x0_y2", "x1_y2", nx_around, True)
add_line("l7", "x1_y2", "x1_y1", ny_foil)
# reuse l2 instead of add_line("l8", "x1_y1", "x0_y1", nx_around, True)

# top left
add_line("l8", "x0_y2", "x0_y3", ny_around, True)
add_line("l9", "x0_y3", "x1_y3", nx_around, True)
add_line("l10", "x1_y3", "x1_y2", ny_around, True)
# reuse l6 instead of re-adding "x1_y2" to "x0_y2"

# bottom mid
# reuse l3: "x1_y0" to "x1_y1"
add_line("l11", "x1_y1", "x2_y1", nx_foil)
add_line("l12", "x2_y1", "x2_y0", ny_around, True)
add_line("l13", "x2_y0", "x1_y0", nx_foil)

# mid mid
add_line("l14", "x1_y2", "x2_y2", nx_foil)
add_line("l15", "x2_y2", "x2_y1", ny_foil)
# reuse l11: "x1_y1" to "x2_y1"
# reuse l7:  "x1_y2" to "x1_y1"

# top mid
add_line("l16", "x1_y3", "x2_y3", nx_foil)
add_line("l17", "x2_y3", "x2_y2", ny_around, True)
# reuse l10: "x1_y3" to "x1_y2"
# reuse l14: "x1_y2" to "x2_y2"

# bottom right
add_line("l18", "x2_y1", "x3_y1", nx_around, True)
add_line("l19", "x3_y1", "x3_y0", ny_around, True)
add_line("l20", "x3_y0", "x2_y0", nx_around, True)
# reuse l12: "x2_y1" to "x2_y0"

# mid right
add_line("l21", "x2_y2", "x3_y2", nx_around, True)
add_line("l22", "x3_y2", "x3_y1", ny_foil)
# reuse l18: "x2_y1" to "x3_y1"
# reuse l15: "x2_y2" to "x2_y1"

# top right
add_line("l23", "x2_y3", "x3_y3", nx_around, True)
add_line("l24", "x3_y3", "x3_y2", ny_around, True)
# reuse l21: "x2_y2" to "x3_y2"
# reuse l17: "x2_y3" to "x2_y2"



s = {}
# Now define 9 transfinite 4-line surfaces
def make_rect(label, linelabels):
    lines = [l[llabel] for llabel in linelabels]
    loop = gmsh.model.geo.addCurveLoop(lines, reorient=True)
    s[label] = gmsh.model.geo.addPlaneSurface([loop])

# left column
make_rect("s1", ["l1", "l2", "l3", "l4"])
make_rect("s2", ["l5", "l6", "l7", "l2"])
make_rect("s3", ["l8", "l9", "l10", "l6"])
# middle column
make_rect("s4", ["l3", "l11", "l12", "l13"])
make_rect("s5", ["l7", "l14", "l15", "l11"])
make_rect("s6", ["l10", "l16", "l17", "l14"])
# right column
make_rect("s7", ["l12", "l18", "l19", "l20"])
make_rect("s8", ["l15", "l21", "l22", "l18"])
make_rect("s9", ["l17", "l23", "l24", "l21"])

# Top: s3 s6 s9
# Mid: s2 s5 s8
# Bot: s1 s4 s7



gmsh.model.geo.synchronize()
def makeTransfiniteCurve(label, start, end, n, prog=0):
    if prog != 0:
        gmsh.model.geo.mesh.setTransfiniteCurve(l[label], n, coef=-prog)
    else:
        gmsh.model.geo.mesh.setTransfiniteCurve(l[label], n)

# bottom left
makeTransfiniteCurve("l1", "x0_y0", "x0_y1", ny_around, prog_y)
makeTransfiniteCurve("l2", "x0_y1", "x1_y1", nx_around, prog_x)
makeTransfiniteCurve("l3", "x1_y1", "x1_y0", ny_around, -prog_y)
makeTransfiniteCurve("l4", "x1_y0", "x0_y0", nx_around, -prog_x)

# mid left
makeTransfiniteCurve("l5", "x0_y1", "x0_y2", ny_foil)
makeTransfiniteCurve("l6", "x0_y2", "x1_y2", nx_around, prog_x)
makeTransfiniteCurve("l7", "x1_y2", "x1_y1", ny_foil)

# top left
makeTransfiniteCurve("l8", "x0_y2", "x0_y3", ny_around, -prog_y)
makeTransfiniteCurve("l9", "x0_y3", "x1_y3", nx_around, prog_x)
makeTransfiniteCurve("l10", "x1_y3", "x1_y2", ny_around, prog_y)

# bottom mid
makeTransfiniteCurve("l11", "x1_y1", "x2_y1", nx_foil)
makeTransfiniteCurve("l12", "x2_y1", "x2_y0", ny_around, -prog_y)
makeTransfiniteCurve("l13", "x2_y0", "x1_y0", nx_foil)

# mid mid
makeTransfiniteCurve("l14", "x1_y2", "x2_y2", nx_foil)
makeTransfiniteCurve("l15", "x2_y2", "x2_y1", ny_foil)

# top mid
makeTransfiniteCurve("l16", "x1_y3", "x2_y3", nx_foil)
makeTransfiniteCurve("l17", "x2_y3", "x2_y2", ny_around, prog_y)

# bottom right
makeTransfiniteCurve("l18", "x2_y1", "x3_y1", nx_around, -prog_x)
makeTransfiniteCurve("l19", "x3_y1", "x3_y0", ny_around, -prog_y)
makeTransfiniteCurve("l20", "x3_y0", "x2_y0", nx_around, prog_x)

# mid right
makeTransfiniteCurve("l21", "x2_y2", "x3_y2", nx_around, -prog_x)
makeTransfiniteCurve("l22", "x3_y2", "x3_y1", ny_foil)

# top right
makeTransfiniteCurve("l23", "x2_y3", "x3_y3", nx_around, -prog_x)
makeTransfiniteCurve("l24", "x3_y3", "x3_y2", ny_around, prog_y)



gmsh.model.geo.synchronize()
for label in s.keys():
    # gmsh.model.mesh.setTransfiniteSurface(s[label], n)
    gmsh.model.geo.mesh.setTransfiniteSurface(s[label])



gmsh.model.geo.synchronize()



# 2. Generate mesh
gmsh.model.geo.synchronize()

inletPointIDs = [p[label] for label in ["x0_y0", "x0_y1", "x0_y2", "x0_y3"]]
gmsh.model.addPhysicalGroup(0, inletPointIDs, tag=300, name="Inlet")
inletLineIDs = [l[label] for label in ["l20", "l13", "l4", "l1", "l5", "l8", "l9", "l16", "l23"]]
print(inletLineIDs)
gmsh.model.addPhysicalGroup(1, inletLineIDs, tag=300, name="Inlet")

outletPointIDs = [p[label] for label in ["x3_y3", "x3_y2", "x3_y1", "x3_y0"]]
gmsh.model.addPhysicalGroup(0, outletPointIDs, tag=302, name="Outlet")
outletLineIDs = [l[label] for label in ["l24", "l22", "l19"]]
gmsh.model.addPhysicalGroup(1, outletLineIDs, tag=302, name="Outlet")

gmsh.model.addPhysicalGroup(2, [s[label] for label in s.keys()], tag=500)



gmsh.model.geo.synchronize()
[gmsh.model.mesh.setRecombine(2, s[label]) for label in s.keys()]
gmsh.model.mesh.generate(2)



# 3. Export to MSH file (temporary)
tmp_msh = "full_mesh.msh"
gmsh.write(tmp_msh)
gmsh.finalize()



# 4. Load with meshio
mesh = meshio.read(tmp_msh)

foil_points = [(x[i], y[i]) for i in range(len(x))]
foil_polygon = shapely.geometry.Polygon(foil_points)

# 5. Filter elements inside the circle
verteces = mesh.cells_dict.get("vertex")
lines = mesh.cells_dict.get("line")
quads = mesh.cells_dict.get("quad")
points = mesh.points

def cell_center(cell):
    coords = points[cell]
    return np.mean(coords[:, :2], axis=0)  # 2D center (x, y)

filtered_entities = []
filtered_quads = []
filtered_lines = []
filtered_verts = []

for vertex in verteces:
    coords = points[vertex]
    center = np.mean(coords[:, :2], axis=0)  # 2D center (x, y)
    p = shapely.geometry.Point(center[0], center[1])
    if not foil_polygon.contains(p):
        filtered_verts.append(meshio.CellBlock("vertex", vertex))

for line in lines:
    coords = points[line]
    center = np.mean(coords[:, :2], axis=0)  # 2D center (x, y)
    p = shapely.geometry.Point(center[0], center[1])
    if not foil_polygon.contains(p):# and domain_polygon.contains(p):
        filtered_lines.append(meshio.CellBlock("line", line))
    # if min(points[cell][0]) <= xmin:

for quad in quads:
    coords = points[quad]
    center = np.mean(coords[:, :2], axis=0)  # 2D center (x, y)
    p = shapely.geometry.Point(center[0], center[1])
    if not foil_polygon.contains(p):
        filtered_quads.append(meshio.CellBlock("quad", quad))

filtered_entities = filtered_quads + filtered_lines + filtered_entities
# 6. overwrite mesh.cells
mesh.cells = filtered_entities
filtered_mesh = meshio.Mesh(
    points=points,
    cells=[("quad", np.array(filtered_quads))],#,("line", np.array(filtered_lines)),("vertex", np.array(filtered_verts))],
    cell_data=None,
    point_data=None,
)

# 7. Save result
filtered_mesh.write("clipped_cartesian_by_type.msh", file_format="gmsh22")
# mesh.write("clipped_cartesian_all.msh", file_format="gmsh22")

# Clean up temp
#   os.remove(tmp_msh)
