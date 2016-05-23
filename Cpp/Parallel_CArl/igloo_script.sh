#!/bin/bash

#PBS -S /bin/bash
#PBS -N carl_test
#PBS -o output_igloo.txt
#PBS -e error_igloo.txt
#PBS -l walltime=00:20:00
#PBS -q iceq
#PBS -l select=1:ncpus=8:mpiprocs=8
#PBS -M thiago.milanetto-schlittler@centralesupelec.fr

# Chargement des modules
module load intel-mpi/5.1.0

# On se place dans le repertoire depuis lequel le job a ete soumis
cd $PBS_O_WORKDIR

mpirun -np 8 ./parallel_intersection_test --inputfile intersection_test_2.inter_test --redirect-stdout --keep-cout
