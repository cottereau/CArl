#PBS -l walltime=0:10:00
#PBS -l select=1:ncpus=4:mpiprocs=4

# "fusion" PBS options
# #PBS -P [PROJECT NAME]

# Charge the modules here
# "fusion" cluster modules
# module purge
# module load intel-compilers/17.0.4
# module load intel-mpi/2017.0.3
# module load intel-mkl/2017.0.3
# module load libgmp/6.1.2
# module load libmpfr/3.1.5

cd $PBS_O_WORKDIR
