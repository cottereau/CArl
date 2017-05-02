#!/bin/bash
echo " -------- Assemble clamped coupled model ..."
mpirun -n $1 ./libmesh_assemble_lin_homogeneous__min_x_clamped -i examples/coupled_traction_test/FETI_solver/assemble_brick_traction_A_1k.txt > assemble_clamped_coupled_model.log

echo " -------- Assemble traction coupled model ..."
mpirun -n $1 ./libmesh_assemble_lin_homogeneous__max_x_traction -i examples/coupled_traction_test/FETI_solver/assemble_brick_traction_B_1k.txt > assemble_traction_coupled_model.log

echo " -------- Assemble coupling matrices ..."
mpirun -n $1 ./CArl_assemble_coupling -i examples/coupled_traction_test/FETI_solver/assemble_coupling_1k.txt > assemble_coupling_matrices.log

echo " -------- FETI setup - initialize ..."
mpirun -n $1 ./CArl_FETI_setup_init -i examples/coupled_traction_test/FETI_solver/setup_FETI_solver_1k.txt

echo " -------- Solve model A - uncoupled ..."
mpirun -n $1 ./libmesh_solve_linear_system -i examples/coupled_traction_test/FETI_solver/scratch_folder/ext_solver_u0_A.txt

echo " -------- Solve model B - uncoupled ..."
mpirun -n $1 ./libmesh_solve_linear_system -i examples/coupled_traction_test/FETI_solver/scratch_folder/ext_solver_u0_B.txt

echo " -------- Solve model A - initialize ..."
mpirun -n $1 ./libmesh_solve_linear_system -i examples/coupled_traction_test/FETI_solver/scratch_folder/ext_solver_A.txt

echo " -------- Solve model B - initialize ..."
mpirun -n $1 ./libmesh_solve_linear_system -i examples/coupled_traction_test/FETI_solver/scratch_folder/ext_solver_B.txt

echo " -------- FETI setup - finalize ..."
mpirun -n $1 ./CArl_FETI_setup_finish -i examples/coupled_traction_test/FETI_solver/scratch_folder/CArl_FETI_setup_finish.txt

echo " -------- Solve model A - iter. no. 0 ..."
mpirun -n $1 ./libmesh_solve_linear_system -i examples/coupled_traction_test/FETI_solver/scratch_folder/ext_solver_A.txt

echo " -------- Solve model B - iter. no. 0 ..."
mpirun -n $1 ./libmesh_solve_linear_system -i examples/coupled_traction_test/FETI_solver/scratch_folder/ext_solver_B.txt