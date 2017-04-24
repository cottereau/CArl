#!/bin/bash

echo " -------- Assemble single model traction ..."
mpirun -n $1 ./libmesh_assemble_lin_homogeneous__traction_test -i examples/single_traction_test/assemble_traction_test_1k.txt > assemble_single_model.log

echo " -------- Solve single model traction ..."
mpirun -n $1 ./libmesh_solve_linear_system -i examples/single_traction_test/solve_traction_test_1k.txt > solve_single_model.log

echo " -------- Assemble clamped coupled model ..."
mpirun -n $1 ./libmesh_assemble_lin_homogeneous__min_x_clamped -i examples/coupled_traction_test/FETI_solver/assemble_brick_traction_A_1k.txt > assemble_clamped_coupled_model.log

echo " -------- Assemble traction coupled model ..."
mpirun -n $1 ./libmesh_assemble_lin_homogeneous__max_x_traction -i examples/coupled_traction_test/FETI_solver/assemble_brick_traction_B_1k.txt > assemble_traction_coupled_model.log

echo " -------- Assemble coupling matrices ..."
mpirun -n $1 ./CArl_assemble_coupling -i examples/coupled_traction_test/FETI_solver/assemble_coupling_1k.txt > assemble_coupling_matrices.log
