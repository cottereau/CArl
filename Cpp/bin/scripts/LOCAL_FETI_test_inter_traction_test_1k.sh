#!/bin/bash

mkdir -p examples/coupled_traction_test/intersection/output/inter_1k
mpirun -np 4 ./CArl_build_intersections -i examples/coupled_traction_test/intersection/inter_traction_test_1k.txt
