#!pvpython
# Program to calculate the energy spectrum of NATriuM output
# Requires:
# 		- NATriuM's fft executable (compiled on your system with fftw3)
#       - paraview (this script has to be executed with Paraview's pvpython)
#       - numpy (for input/output of tables)

import sys
import os
import numpy
#from vtk import *
from paraview.simple import *


print "Usage: pvpython compute_energy_spectrum.py <vtk_file>"

argc = len (sys.argv)
vtu_file = sys.argv[1]
print "VTU file: ", vtu_file

# Read file
vtu_reader = OpenDataFile(vtu_file)
vtu_reader.UpdatePipeline()
print vtu_reader.PointData[:]

# Definition file for structured grid
print "Transfer data to structured grid..."
vtk_filename = "foo.vtk"
vtk_file = open(vtk_filename,"w+")
vtk_file.write("# vtk DataFile Version 3.0\n")
vtk_file.write("vtk output\n")
vtk_file.write("ASCII\n")
vtk_file.write("DATASET STRUCTURED_POINTS\n")
vtk_file.write("DIMENSIONS 1000 1000 1\n")
vtk_file.write("SPACING 0.001 0.001 0.001\n")
vtk_file.write("ORIGIN 0.0 0.0 0.00\n\n")
vtk_file.close()
vtk_reader = OpenDataFile(vtk_filename)
vtk_reader.UpdatePipeline()


# Interpolate on structured grid
filter = ResampleWithDataset(Source=vtk_reader, Input=vtu_reader)
filter.UpdatePipeline()
print "done."

# Write to CSV
writer = CreateWriter("foo.csv", filter)
writer.FieldAssociation = "Points"
writer.UpdatePipeline()
del writer


# Open CSV
