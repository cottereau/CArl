#!/bin/bash

# - Build the intersections
. scripts/LOCAL_Time_Mono_test_inter_traction_dyn_test_Barre.sh
sleep 5

# - Prepare the external solver
. scripts/LOCAL_Time_Mono_test_init_ext_solver_traction_dyn_test_Barre.sh
sleep 5

# - Build the coupling matrices
. scripts/LOCAL_Time_Mono_test_assemble_coupling_matrices_traction_dyn_test_Barre.sh