#!/bin/bash

mkdir -p ../examples/coupled_traction_dyn_test/intersection/output/inter_Barre
mpirun -np 1 ./CArl_build_intersections -i ../examples/coupled_traction_dyn_test/intersection/inter_traction_test_Barre.txt
