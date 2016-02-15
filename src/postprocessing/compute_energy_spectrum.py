#!pvpython
# -*- coding: utf-8 -*-
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
from math import *

# dimension of flow domain
dim = 2
# spacing between structured grid points
spacing = 0.1
# number of grid points in each direction (x,y,z)
N = [101, 31, 1]
# origin of structured grid discretization
origin = [10.0, 1.0, 0.0]
TEST = False 
#volume of region for which the fft is calculated
volume_element = pow( spacing, 3)
PI = 4*atan(1)

if not TEST:
    print "Usage: pvpython compute_energy_spectrum.py <vtu_file>"
    
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
    vtk_file.write("DIMENSIONS %i %i %i\n" %(N[0], N[1], N[2]) )
    vtk_file.write("SPACING %f %f %f\n" %(spacing,spacing,spacing))
    vtk_file.write("ORIGIN %f %f %f\n\n" %(origin[0], origin[1], origin[2]) )
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
    
    # Split CSV files into velocity components
    A = numpy.loadtxt("foo.csv", delimiter=',', skiprows=1)
    col_ux = 1
    col_uy = 2
    col_uz = 3
    numpy.savetxt( "ux.txt", A[:,col_ux])
    numpy.savetxt( "uy.txt", A[:,col_uy])
    if dim == 3:
        numpy.savetxt( "uz.txt", A[:,col_uz])

if TEST:
    test_array = []
    for i in range(N[0]):
        for j in range(N[1]):
            for l in range(N[2]):
                test_array += [[cos(8*atan(1)*(i)/N[0])  + sin(8*atan(1)*(4.0*j)/N[1]), cos(8*atan(1)*(2.0*j)/N[1]) + sin(8*atan(1)*(12.0 *j)/N[1]) , 0.0]]
    test_array = numpy.array(test_array)
    numpy.savetxt( "ux.txt", test_array[:,0])
    numpy.savetxt( "uy.txt", test_array[:,1])
    numpy.savetxt( "uz.txt", test_array[:,2])

# Fourier transform
print "Fourier transform of velocity components..."
script = os.path.join(os.getenv("NATRIUM_DIR"),"src/postprocessing/fft")
os.system("%s %d %d %d %d %s %s" %(script, dim, N[0], N[1], N[2], "ux.txt", "ux_hat.txt"))
os.system("%s %d %d %d %d %s %s" %(script, dim, N[0], N[1], N[2], "uy.txt", "uy_hat.txt"))
if (dim == 3):
    os.system("%s %d %d %d %d %s %s" %(script, dim, N[0], N[1], N[2], "uz.txt", "uz_hat.txt"))
print "done."

# Load fourier coefficients
ux_hat = numpy.loadtxt("ux_hat.txt",delimiter=' ',  skiprows=1)
kx  = ux_hat[:,0]
ky  = ux_hat[:,1]
kz  = ux_hat[:,2]
ux_hat = ux_hat[:,3:]
uy_hat = numpy.loadtxt("uy_hat.txt",delimiter=' ',  skiprows=1)[:,3:]
if (dim == 3):
    uz_hat = numpy.loadtxt("uz_hat.txt",delimiter=' ',  skiprows=1)[:,3:]
# (two columns: Re | Im)


# Assemble energy spectrum
n_tics  = int( max(N))
sum_weight = [0.0 for i in range(n_tics)]
energies = [0.0 for i in range(n_tics)]
sigma = 0.1
max_kappa = sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2])

for i in range(N[0]):
    for j in range(N[1]):
        for l in range(N[2]):
            # calculate wave numbers and kappa
            # center fourier modes around 0
            ik = 0.5 * N[0] - fabs( i - 0.5*N[0])
            jk = 0.5 * N[1] - fabs( j - 0.5*N[1])
            lk = 0.5 * N[2] - fabs( l - 0.5*N[2])
            kappa = sqrt(ik**2 + jk**2 + lk**2)
            pos_k = N[1]*N[2] * i + N[2] * j + l
            assert kx[pos_k] == i
            assert ky[pos_k] == j
            assert kz[pos_k] == l
            E = 0.5 * (ux_hat[pos_k][0] * ux_hat[pos_k][0] + ux_hat[pos_k][1] * ux_hat[pos_k][1] + uy_hat[pos_k][0] * uy_hat[pos_k][0] + uy_hat[pos_k][1] * uy_hat[pos_k][1] )
            if dim == 3:
                E +=  0.5* ( uz_hat[pos_k][0] * uz_hat[pos_k][0] + uz_hat[pos_k][1] * uz_hat[pos_k][1] )
            k_index = 1.0 * kappa / max_kappa * n_tics
            # add to neighoring modes
            for m in range( int (max(0,kappa-5)) , int( min(n_tics,kappa+5) ) ) :
                sum_weight[m] = sum_weight[m] + exp( -0.5 * ((1.0*k_index - m)/sigma)**2 )
                energies[m] = energies[m] + exp( -0.5 *((1.0*k_index - m)/sigma)**2 ) * E

# calculate real wave numbers
#wave_numbers  = [ 2.0*PI*i /   ]
wave_numbers = [i for i in range(n_tics)]
for i in range(n_tics):
    if (sum_weight[i] > 1e-40):
        energies[i] /= sum_weight[i]
        # from discrete FT to FT
        energies[i] *= (volume_element * pow(sqrt(2 * PI), dim) )


# rescale to physical quantities
#real_phases

numpy.savetxt("E.txt", [[wave_numbers[i], energies[i]] for i in range(n_tics)])
print "Calculation finished. Energy spectrum in E.txt"
