#!/bin/bash
# NOTE : remember to set the modules and any other PBS options before launching

# - Build the intersections
job_inter=`qsub scripts/PBS_FETI_test_inter_traction_test_8k.pbs`

# - Prepare the external solver
job_ext=`qsub scripts/PBS_FETI_test_assemble_ext_solver_traction_test_8k.pbs`

# - Build the coupling matrices AFTER finishing $job_inter
job_coupl=`qsub -W depend=afterok:$job_inter scripts/PBS_FETI_test_assemble_coupling_matrices_traction_test_8k.pbs`

# - Solve the coupled system AFTER finishing both $job_ext and $job_coupl
job_solve=`qsub -W depend=afterok:$job_ext:$job_coupl scripts/PBS_FETI_launch_coupled_solver_8k.pbs`

# - Set the coupled solution mesh AFTER finishing $job_solve
echo " ---> After all the solver jobs are finished, execute:"
echo "      qsub scripts/PBS_FETI_apply_solution_traction_test_8k.pbs"
