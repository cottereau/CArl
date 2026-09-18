#!/bin/bash
# NOTE : remember to set the modules and any other PBS options before launching

# - Build the intersections
job_inter=`qsub ../examples/scripts/PBS_Time_Mono_test_inter_traction_dyn_test.pbs`

# - Prepare the external solver
job_ext=`qsub ../examples/scripts/PBS_Time_Mono_test_init_ext_solver_traction_dyn_test.pbs`

# - Build the coupling matrices AFTER finishing $job_inter
job_coupl=`qsub -W depend=afterok:$job_inter ../examples/scripts/PBS_Time_Mono_test_assemble_coupling_matrices_traction_dyn_test.pbs`

# - Solve the coupled system AFTER finishing both $job_ext and $job_coupl
job_solve=`qsub -W depend=afterok:$job_ext:$job_coupl ../examples/scripts/PBS_Time_Mono_FETI_launch_coupled_solver.pbs`
