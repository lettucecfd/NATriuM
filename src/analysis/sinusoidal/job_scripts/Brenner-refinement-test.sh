# the following special lines starting with #PBS tells the batch system what resources you need
##########PBS -N Brenner-Ma
#PBS -q default
#PBS -l nodes=1:ppn=1
#PBS -l walltime=12:00:00
#PBS -l vmem=8GB
#PBS -V

module load gcc/4.9.0
module load openmpi/gnu

cp /home/akraem3m/NATriuM/job_scripts/compare-Brenner/refinement-test-2015-06-19/plot_flowfactors_generic.gnp $OUTPUT_DIR
cd $OUTPUT_DIR
echo "CFG_ID="$CFG_ID", Ma ="$Ma", Gamma="$Gamma", Refinement="$Refinement",orderFE="$orderFE",cfl="$CFL
time /home/akraem3m/NATriuM/bin_release/examples/step-sinusoidal/Brenner-refinement-test $CFG_ID $Ma $Gamma $Refinement $orderFE $CFL
gnuplot plot_flowfactors_generic.gnp
echo "Plots created. Done."

