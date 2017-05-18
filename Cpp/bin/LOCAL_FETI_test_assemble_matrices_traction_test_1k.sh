#!/bin/bash
mkdir -p examples/coupled_traction_test/FETI_solver/brick_traction_1k/system_matrices

echo " --- Assemble clamped coupled model ..."
mpirun -n 4 ./libmesh_assemble_lin_homogeneous__min_x_clamped -i examples/coupled_traction_test/FETI_solver/brick_traction_1k/assemble_brick_traction_A_1k.txt > assemble_clamped_coupled_model.log

echo " --- Assemble traction coupled model ..."
mpirun -n 4 ./libmesh_assemble_lin_homogeneous__max_x_traction -i examples/coupled_traction_test/FETI_solver/brick_traction_1k/assemble_brick_traction_B_1k.txt > assemble_traction_coupled_model.log

echo " --- Assemble coupling matrices ..."
mpirun -n 4 ./CArl_assemble_coupling -i examples/coupled_traction_test/FETI_solver/brick_traction_1k/assemble_coupling_1k.txt > assemble_coupling_matrices.log
