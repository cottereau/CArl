#!/bin/bash


#SBATCH --job-name=test_s_1k
#SBATCH --time=03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#
# "ruche" PBS options


# Charge the modules here
# "ruche" cluster modules
 #module purge
 #module load intel-compilers/17.0.4
 #module load intel-mpi/2017.0.3
 #module load intel-mkl/2017.0.3
 #module load libgmp/6.1.2
 #module load libmpfr/3.1.5

# On se place dans le repertoire depuis lequel le job a ete soumis
cd ${SLURM_SUBMIT_DIR}

srun -n 4 ./libmesh_apply_solution_homogeneous -i ./FETI_solver/brick_traction_1k/apply_solution_brick_traction_A_1k.txt

srun -n 4 ./libmesh_apply_solution_homogeneous -i ./FETI_solver/brick_traction_1k/apply_solution_brick_traction_B_1k.txt  
