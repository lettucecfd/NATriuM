import numpy
import math

A = numpy.loadtxt("../results_table.txt")
B = numpy.loadtxt("../global_turbulence_table.txt")
derivatives = []

out_table = [] 

n = 5

for i in range (len(A)-(n+1)):
  #compressible
  derivative = (A[i+n][4] - A[i][4]) /(A[i+n][1] - A[i][1]) 
  normalized = derivative / ((8*math.atan(1))**3) #/m^3
  #incompressible
  kE_i = 0.5 * (B[i][32]+B[i][38]+B[i][45])
  kE_i_plus_n = 0.5 * (B[i+n][32]+B[i+n][38]+B[i+n][45])
  derivative2 = (kE_i_plus_n-kE_i)/(A[i+n][1] - A[i][1]) 
  normalized2 = derivative2 / ((8*math.atan(1))**3)
  out_table = out_table + [[A[i][1], normalized, normalized2]]

C = numpy.array(out_table)
numpy.savetxt("derivatives.txt", C)


