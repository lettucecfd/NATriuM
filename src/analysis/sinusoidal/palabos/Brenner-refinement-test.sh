# the following special lines starting with #PBS tells the batch system what resources you need
##PBS -N AbruchTest 700
#PBS -q default
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -l vmem=2GB
#PBS -V

module load gcc/4.9.0
module load openmpi/gnu

cd /home/akraem3m/bin/palabos-v1.5r1/vincent_tobias
# Ma = 0.01
time mpirun -np 1 ./projektwoche1 $Refinement 0.005773503 $ChannelLength $ChannelAmplitude 3.333333

