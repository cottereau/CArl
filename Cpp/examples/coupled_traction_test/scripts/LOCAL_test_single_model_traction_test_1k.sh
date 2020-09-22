#!/bin/bash

echo " --- Assemble single model traction ..."
mpirun -n 4 ./libmesh_assemble_lin_homogeneous__traction_test -i $CARLBUILD/examples/single_traction_test/assemble_traction_test_1k.txt > assemble_single_model.log

echo " --- Solve single model traction ..."
mpirun -n 4 ./libmesh_solve_linear_system -i $CARLBUILD/examples/single_traction_test/solve_traction_test_1k.txt > solve_single_model.log
