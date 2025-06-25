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
dx_foil = 0.05
dx_around = 0.2

# foil coordinates
x = [1.0, 0.999748, 0.998993, 0.997736, 0.995977, 0.993719, 0.990964, 0.987715, 0.983974, 0.979746, 0.975036, 0.969846, 0.964184, 0.958054, 0.951463, 0.944418, 0.936925, 0.928992, 0.920627, 0.911838, 0.902635, 0.893027, 0.883022, 0.872632, 0.861867, 0.850737, 0.839255, 0.82743, 0.815276, 0.802805, 0.790028, 0.77696, 0.763613, 0.75, 0.736136, 0.722033, 0.707708, 0.693173, 0.678443, 0.663534, 0.64846, 0.633237, 0.617879, 0.602403, 0.586824, 0.571157, 0.555419, 0.539625, 0.523791, 0.507933, 0.492067, 0.476209, 0.460375, 0.444581, 0.428843, 0.413176, 0.397597, 0.382121, 0.366763, 0.35154, 0.336466, 0.321557, 0.292292, 0.277967, 0.263864, 0.25, 0.236387, 0.22304, 0.209972, 0.197195, 0.184724, 0.17257, 0.160745, 0.149263, 0.138133, 0.127368, 0.116978, 0.106973, 0.0973649, 0.0881617, 0.0793732, 0.0710083, 0.0630753, 0.0555823, 0.0485367, 0.0419458, 0.035816, 0.0301537, 0.0249644, 0.0202535, 0.0160256, 0.0122851, 0.00903565, 0.00628056, 0.00402259, 0.00226404, 0.00100666, 0.000251729, 0.0, 0.000251729, 0.00100666, 0.00226404, 0.00402259, 0.00628056, 0.00903565, 0.0122851, 0.0160256, 0.0202535, 0.0249644, 0.0301537, 0.035816, 0.0419458, 0.0485367, 0.0555823, 0.0630753, 0.0710083, 0.0793732, 0.0881617, 0.0973649, 0.106973, 0.116978, 0.127368, 0.138133, 0.149263, 0.160745, 0.17257, 0.184724, 0.197195, 0.209972, 0.22304, 0.236387, 0.25, 0.263864, 0.277967, 0.306827, 0.321557, 0.336466, 0.35154, 0.366763, 0.382121, 0.397597, 0.413176, 0.428843, 0.444581, 0.460375, 0.476209, 0.492067, 0.507933, 0.523791, 0.539625, 0.555419, 0.571157, 0.586824, 0.602403, 0.617879, 0.633237, 0.64846, 0.663534, 0.678443, 0.693173, 0.707708, 0.722033, 0.736136, 0.75, 0.763613, 0.77696, 0.790028, 0.802805, 0.815276, 0.82743, 0.839255, 0.850737, 0.861867, 0.872632, 0.883022, 0.893027, 0.902635, 0.911838, 0.920627, 0.928992, 0.936925, 0.944418, 0.951463, 0.958054, 0.964184, 0.969846, 0.975036, 0.979746, 0.983974, 0.987715, 0.990964, 0.993719, 0.995977, 0.997736, 0.998993, 0.999748];
y = [-1.66533e-17, 3.65828e-05, 0.000146223, 0.000328595, 0.00058316, 0.00090917, 0.00130567, 0.00177151, 0.00230534, 0.00290565, 0.00357073, 0.00429874, 0.00508767, 0.00593537, 0.00683958, 0.00779793, 0.00880793, 0.00986703, 0.0109726, 0.0121219, 0.0133121, 0.0145406, 0.0158044, 0.0171006, 0.0184265, 0.0197789, 0.0211551, 0.0225521, 0.0239669, 0.0253966, 0.0268382, 0.0282887, 0.0297451, 0.0312044, 0.0326635, 0.0341194, 0.035569, 0.0370091, 0.0384366, 0.0398482, 0.0412407, 0.0426107, 0.0439549, 0.0452698, 0.0465521, 0.0477982, 0.0490045, 0.0501675, 0.0512835, 0.0523489, 0.0533601, 0.0543134, 0.0552053, 0.0560321, 0.0567903, 0.0574765, 0.0580873, 0.0586193, 0.0590696, 0.059435, 0.0597128, 0.0599003, 0.0599951, 0.0598982, 0.0597028, 0.0594075, 0.0590111, 0.0585128, 0.0579121, 0.0572087, 0.0564028, 0.0554948, 0.0544854, 0.0533756, 0.0521668, 0.0508605, 0.0494586, 0.0479633, 0.0463768, 0.0447018, 0.042941, 0.0410972, 0.0391735, 0.037173, 0.0350989, 0.0329544, 0.0307426, 0.0284669, 0.0261302, 0.0237357, 0.0212862, 0.0187844, 0.0162331, 0.0136345, 0.0109908, 0.008304, 0.0055757, 0.00280732, 0.0, -0.00280732, -0.0055757, -0.008304, -0.0109908, -0.0136345, -0.0162331, -0.0187844, -0.0212862, -0.0237357, -0.0261302, -0.0284669, -0.0307426, -0.0329544, -0.0350989, -0.037173, -0.0391735, -0.0410972, -0.042941, -0.0447018, -0.0463768, -0.0479633, -0.0494586, -0.0508605, -0.0521668, -0.0533756, -0.0544854, -0.0554948, -0.0564028, -0.0572087, -0.0579121, -0.0585128, -0.0590111, -0.0594075, -0.0597028, -0.0598982, -0.0599951, -0.0599003, -0.0597128, -0.059435, -0.0590696, -0.0586193, -0.0580873, -0.0574765, -0.0567903, -0.0560321, -0.0552053, -0.0543134, -0.0533601, -0.0523489, -0.0512835, -0.0501675, -0.0490045, -0.0477982, -0.0465521, -0.0452698, -0.0439549, -0.0426107, -0.0412407, -0.0398482, -0.0384366, -0.0370091, -0.035569, -0.0341194, -0.0326635, -0.0312044, -0.0297451, -0.0282887, -0.0268382, -0.0253966, -0.0239669, -0.0225521, -0.0211551, -0.0197789, -0.0184265, -0.0171006, -0.0158044, -0.0145406, -0.0133121, -0.0121219, -0.0109726, -0.00986703, -0.00880793, -0.00779793, -0.00683958, -0.00593537, -0.00508767, -0.00429874, -0.00357073, -0.00290565, -0.00230534, -0.00177151, -0.00130567, -0.00090917, -0.00058316, -0.000328595, -0.000146223, -3.65828e-05];

foil_top = max(y)
foil_bottom = min(y)
foil_front = min(x)
foil_back = max(x)

# 1. Define rectangular Cartesian mesh
inlet_top = gmsh.model.geo.addPoint(xmin, ymax, 0)
inlet_foil_top = gmsh.model.geo.addPoint(xmin, foil_back, 0)
inlet_foil_bottom = gmsh.model.geo.addPoint(xmin, foil_front, 0)
inlet_bottom = gmsh.model.geo.addPoint(xmin, ymin, 0)

foil_front_bottom = gmsh.model.geo.addPoint(foil_bottom, ymin, 0)
foil_back_bottom = gmsh.model.geo.addPoint(foil_top, ymin, 0)

outlet_bottom = gmsh.model.geo.addPoint(xmax, ymax, 0)
outlet_foil_bottom = gmsh.model.geo.addPoint(xmax, foil_front, 0)
outlet_foil_top = gmsh.model.geo.addPoint(xmax, foil_back, 0)
outlet_top = gmsh.model.geo.addPoint(xmax, ymax, 0)

foil_back_top = gmsh.model.geo.addPoint(foil_top, ymax, 0)
foil_front_top = gmsh.model.geo.addPoint(foil_bottom, ymax, 0)

domain_points = [inlet_top, inlet_foil_top, inlet_foil_bottom, inlet_bottom,
                foil_front_bottom, foil_back_bottom, outlet_bottom, outlet_foil_bottom, outlet_foil_top, outlet_top,
                foil_back_top, foil_front_top]
# [inlet_top 0 inlet_foil_top 1 inlet_foil_bottom 2 inlet_bottom,
# 3 foil_front_bottom 4 foil_back_bottom 5 outlet_bottom 6 outlet_foil_bottom 7 outlet_foil_top 8 outlet_top
# 9 foil_back_top 10 foil_front_top 11]

# Top row:    [s1] [s2] [s3]
# Middle row: [s4] [s5] [s6]
# Bottom row: [s7] [s8] [s9]

lines = []
for i in range(len(domain_points)-1):
    lines.append(gmsh.model.geo.addLine(domain_points[i], domain_points[i+1]))
lines.append(gmsh.model.geo.addLine(domain_points[-1], domain_points[0]))

domain_loop = gmsh.model.geo.addCurveLoop(lines)
domain_surface = gmsh.model.geo.addPlaneSurface([domain_loop])

gmsh.model.geo.synchronize()

# Transfinite (structured) meshing
nx_foil = int((max(x)-min(x)) / dx_foil) + 1
ny_foil = int((max(y)-min(y)) / dx_foil) + 1
nx_around = int(Lx / dx_around) + 1
ny_around = int(Lx / dx_around) + 1

print(lines)

for i in [0, 2, 6, 8]:
    gmsh.model.mesh.setTransfiniteCurve(lines[i], nx_around, meshType="Progression", coef=1.1)
for i in [3, 5, 9, 11]:
    gmsh.model.mesh.setTransfiniteCurve(lines[i], ny_around, meshType="Progression", coef=1.1)
for i in [1, 7]:
    gmsh.model.mesh.setTransfiniteCurve(lines[i], nx_foil)
for i in [4, 10]:
    gmsh.model.mesh.setTransfiniteCurve(lines[i], ny_foil)
gmsh.model.mesh.setTransfiniteSurface(domain_surface)
gmsh.model.mesh.setRecombine(2, domain_surface)  # make quads

# 2. Generate mesh
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
cells = mesh.cells_dict.get("quad", [])  # or "triangle" if no Recombine
points = mesh.points

def cell_center(cell):
    coords = points[cell]
    return np.mean(coords[:, :2], axis=0)  # 2D center (x, y)

filtered_cells = []
for cell in cells:
    center = cell_center(cell)
    p = shapely.geometry.Point(center[0], center[1])
    if not foil_polygon.contains(p):# and domain_polygon.contains(p):
        filtered_cells.append(cell)

# 6. Create filtered mesh
filtered_mesh = meshio.Mesh(
    points=points,
    cells=[("quad", np.array(filtered_cells))],
    cell_data=None,
    point_data=None,
)

# 7. Save result
filtered_mesh.write("clipped_cartesian.msh", file_format="gmsh22")

# Clean up temp
os.remove(tmp_msh)
