#PBS -l walltime=0:10:00
#PBS -q haswellq
#PBS -P couest
#PBS -l select=1:ncpus=4:mpiprocs=4

# Charge the modules here
# "fusion" cluster modules
# module purge
# module load intel-compilers/16.0.3
# module load intel-mpi/5.1.2

cd $PBS_O_WORKDIR
