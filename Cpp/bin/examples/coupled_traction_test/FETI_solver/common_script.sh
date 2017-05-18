#PBS -l walltime=0:10:00
#PBS -q haswellq
#PBS -P couest
#PBS -l select=1:ncpus=4:mpiprocs=4

# Charge the modules here
module purge
module load boost/1.61.0
module load gcc/5.4.0
module load intel-compilers/16.0.3
module load intel-mkl/11.3.3
module load intel-mpi/5.1.2
module load metis/5.1.0
module load parmetis/4.0.3
module load hdf5/1.8.16

cd $PBS_O_WORKDIR
