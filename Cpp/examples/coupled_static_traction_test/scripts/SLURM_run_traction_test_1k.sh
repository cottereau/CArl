#!/bin/bash
# NOTE : remember to set the modules and any other PBS options before launching

# - Build the intersections
job_inter=`sbatch ./run/SLURM_FETI_test_inter_traction_test_1k.slurm`

# - Prepare the external solver (libMesh)
job_ext=`sbatch ./run/SLURM_FETI_test_assemble_ext_solver_traction_test_1k.slurm`

# - Build the coupling matrices AFTER finishing $job_inter
job_coupl=`sbatch --dependency=afterok:$job_inter ./run/SLURM_FETI_test_assemble_coupling_matrices_traction_test_1k.slurm`

# - Solve the coupled system AFTER finishing both $job_ext and $job_coupl
job_solve=`sbatch --dependency=afterok:$job_ext:$job_coupl ./run/SLURM_FETI_launch_coupled_solver.slurm`

#Comment
# - Set the coupled solution mesh AFTER finishing $job_solve
echo " ---> After all the solver jobs are finished, execute:"
echo "      sbatch ./run/SLURM_FETI_apply_solution_traction_test_1k.slurm"
