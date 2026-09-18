#!/bin/bash
echo " --- Assemble coupling matrices ..."
mpirun --allow-run-as-root -n 1 ./CArl_assemble_coupling -i ../examples/coupled_traction_test/FETI_solver/brick_traction_1k/assemble_coupling_1k.txt > assemble_coupling_matrices.log
