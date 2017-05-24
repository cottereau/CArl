#!/bin/bash

mpirun -n 4 ./libmesh_apply_solution_homogeneous -i examples/coupled_traction_test/FETI_solver/brick_traction_1k/apply_solution_brick_traction_A_1k.txt

mpirun -n 4 ./libmesh_apply_solution_homogeneous -i examples/coupled_traction_test/FETI_solver/brick_traction_1k/apply_solution_brick_traction_B_1k.txt  
