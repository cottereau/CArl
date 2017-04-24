#PBS -l walltime=3:59:00
#PBS -q haswellq
#PBS -P couest
#PBS -l select=16:ncpus=24:mpiprocs=24
#PBS -m abe

# Chargement des modules
module purge
module load boost/1.61.0
module load gcc/5.4.0
module load intel-compilers/16.0.3
module load intel-mkl/11.3.3
module load intel-mpi/5.1.2
module load metis/5.1.0
module load parmetis/4.0.3
module load hdf5/1.8.16

# On se place dans le repertoire depuis lequel le job a ete soumis
cd $PBS_O_WORKDIR