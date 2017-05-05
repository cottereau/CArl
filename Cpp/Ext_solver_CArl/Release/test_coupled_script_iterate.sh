#!/bin/bash
echo " -------- FETI iteration ..."
mpirun -n $1 ./CArl_FETI_iterate -i examples/coupled_traction_test/FETI_solver/scratch_folder/CArl_FETI_iterate.txt

echo " -------- Solve model A - iter. no. 0 ..."
mpirun -n $1 ./libmesh_solve_linear_system -i examples/coupled_traction_test/FETI_solver/scratch_folder/ext_solver_A.txt

echo " -------- Solve model B - iter. no. 0 ..."
mpirun -n $1 ./libmesh_solve_linear_system -i examples/coupled_traction_test/FETI_solver/scratch_folder/ext_solver_B.txt