#!/bin/bash
echo " --- Creating directories ..." 
mkdir -p ../examples/coupled_traction_dyn_test/Time_Mono_solver/Barre_traction/system_matrices
mkdir -p ../examples/coupled_traction_dyn_test/Time_Mono_solver/Barre_traction/system_initial_conditions

echo " --- Initializing external solver for clamped coupled model ..." 
mpirun -n 1 ./libmesh_init_dyn_lin_homogeneous__min_x_clamped -i ../examples/coupled_traction_dyn_test/Time_Mono_solver/Barre_traction/init_barre_traction_A.txt > init_clamped_coupled_model.log

echo " --- Initializing external solver for traction coupled model ..."
mpirun -n 1 ./libmesh_init_dyn_lin_homogeneous__max_x_traction -i ../examples/coupled_traction_dyn_test/Time_Mono_solver/Barre_traction/init_barre_traction_B.txt > init_traction_coupled_model.log
