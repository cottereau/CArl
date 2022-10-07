#!/bin/bash

#SBATCH --job-name=test_si_1k
#SBATCH --time=03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=4


# Charge the modules here
# "fusion" cluster modules
# module purge
# module load intel-compilers/17.0.4
# module load intel-mpi/2017.0.3
# module load intel-mkl/2017.0.3
# module load libgmp/6.1.2
# module load libmpfr/3.1.5

# On se place dans le repertoire depuis lequel le job a ete soumis
cd ${SLURM_SUBMIT_DIR}

srun -np 4 $CARLBUILD/CArl_FETI_setup_init -i ./FETI_solver/brick_traction_1k/PBS_setup_FETI_solver_1k.txt
