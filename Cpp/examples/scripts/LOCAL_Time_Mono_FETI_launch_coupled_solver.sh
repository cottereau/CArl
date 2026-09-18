#!/bin/bash

#mpirun -np -n 4 ./CArl_FETI_setup_init -i ../examples/coupled_traction_test/FETI_solver/brick_traction_1k/LOCAL_setup_FETI_solver_1k.txt
mpirun -n 4 ./CArl_Time_Mono_setup -i ../examples/coupled_traction_dyn_test/Time_Mono_solver/Barre_traction/LOCAL_setup_Time_Mono_solver.txt
