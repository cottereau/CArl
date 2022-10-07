#!/bin/bash

# - Build the intersections
. scripts/LOCAL_FETI_test_inter_traction_test_1k.sh
sleep 5

# - Prepare the external solver
. scripts/LOCAL_FETI_test_assemble_ext_solver_traction_test_1k.sh
sleep 5

# - Build the coupling matrices
. scripts/LOCAL_FETI_test_assemble_coupling_matrices_traction_test_1k.sh
