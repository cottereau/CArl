#!/bin/bash
echo " --- Assemble coupling matrices ..."
mpirun -n 1 ./CArl_assemble_coupling -i ../examples/coupled_traction_dyn_test/Time_Mono_solver/Barre_traction/assemble_coupling_barre.txt > assemble_coupling_matrices.log
