
#!/bin/bash
# NOTE : remember to set the modules and any other PBS options before launching

# - Build the intersections
#job_inter=`qsub ../examples/coupled_dynamic/scripts/PBS_FETI_test_inter_dynamic_1k.pbs`

# - Prepare the external solver
job_ext=`qsub ../examples/coupled_dynamic/scripts/PBS_FETI_test_assemble_ext_solver_traction_test_1k_dyn.pbs`
#
## - Build the coupling matrices AFTER finishing $job_inter
#job_coupl=`qsub -W depend=afterok:$job_inter ../examples/coupled_traction_test/scripts/PBS_FETI_test_assemble_coupling_matrices_traction_test_1k.pbs`
#
## - Solve the coupled system AFTER finishing both $job_ext and $job_coupl
#job_solve=`qsub -W depend=afterok:$job_ext:$job_coupl ../examples/coupled_traction_test/scripts/PBS_FETI_launch_coupled_solver.pbs`
#
## - Set the coupled solution mesh AFTER finishing $job_solve
#echo " ---> After all the solver jobs are finished, execute:"
#echo "      qsub ../examples/coupled_traction_test/scripts/PBS_FETI_apply_solution_traction_test_1k.pbs"
