#!/bin/bash

mkdir -p ../examples/coupled_traction_test/intersection/output/inter_1k
mpirun --allow-run-as-root -np 1 ./CArl_build_intersections -i ../examples/coupled_traction_test/intersection/inter_traction_test_1k.txt
