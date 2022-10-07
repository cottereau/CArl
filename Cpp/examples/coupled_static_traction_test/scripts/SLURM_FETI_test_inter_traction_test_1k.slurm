#!/bin/bash

#SBATCH --job-name=test_i_1ki
#SBATCH --output=ro.txt
#SBATCH --error=re.txt
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

mkdir -p ./intersection/output/inter_1k

srun -np 4 $CARLBUILD/CArl_build_intersections -i ./intersection/inter_traction_test_1k.txt
