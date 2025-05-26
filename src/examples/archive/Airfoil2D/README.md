# README

This example requires some previous steps:

0. Move to `mesh` directory.
1. Create the airfoil point in gmsh. I used https://github.com/jte0419/NACA_4_Digit_Airfoil with 100 
points and open trailing edge (otherwise, the last line will be too small)
2. Use naca2domain.py to create the domain:
`python naca2domain.py`
3. Obtain gmsh from https://gmsh.info/#Download
4. Create the mesh:
```
for angle in 0 1 2 3 4 5
do
    ~/Downloads/gmsh-4.11.1-Linux64/bin/gmsh nonUni_closed/NACA0012_${angle}deg.geo -save
done
```
5. Mind the naming of mesh files (see DiamondObstacle2D::makeGrid)