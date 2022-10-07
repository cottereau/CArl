#!/bin/bash

#SBATCH --job-name=test_ae_1k
#SBATCH --time=03:59:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --output=ro.txt
#SBATCH --error=re.txt


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

mkdir -p ./FETI_solver/brick_traction_1k/system_matrices

echo " --- Assemble clamped coupled model ..."
srun -np 4 $CARLBUILD/libmesh_assemble_lin_homogeneous__min_x_clamped -i ./FETI_solver/brick_traction_1k/assemble_brick_traction_A_1k.txt 

echo " --- Assemble traction coupled model ..."
srun -np 4 $CARLBUILD/libmesh_assemble_lin_homogeneous__max_x_traction -i ./FETI_solver/brick_traction_1k/assemble_brick_traction_B_1k.txt

