#!/bin/bash
echo " --- Assemble coupling matrices ..."
mpirun -n 4 ./CArl_assemble_coupling -i ../examples/coupled_traction_test/FETI_solver/brick_traction_1k/assemble_coupling_1k.txt > assemble_coupling_matrices.log
