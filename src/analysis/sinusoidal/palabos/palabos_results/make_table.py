import os
import numpy
import matplotlib

#create tables from single files
os.system("grep user *cfg?.e* > times.txt")
os.system("grep Flowfactor: *.o* > results.txt")
#load both files and bring into right format
t = numpy.loadtxt("times.txt", usecols=(1,), dtype="S7", unpack=True)
times = []
for i in range(len(t)):
  m,s = t[i].split('m')
  times += [60*float(m)+float(s)]
f = open("results.txt")
res = f.readlines()
table = []
for i in range(len(res)):
  a = res[i].replace('Brenner-Ref','')
  a = a.replace('-cfg',' ')
  a = a.replace('.o', ' ')
  a = a.replace(':Flowfactor:', '')
  ref,cfg,jobid,flowfactor = a.split()
  table += [[ int(ref), int(cfg), float(flowfactor), times[i]]]
#sort
table=sorted(table,key=lambda x:x[0])
#save file
numpy.savetxt("table.txt", table, fmt=['%i','%i','%1.5e','%1.3e'])#, header='# N config Psi_s time(user)')
