#!/bin/bash
# Run example named "test" on 2 nodes - 24 MPI processes (12 MPI processes per node):
#
#   ./prepare_batch_run.sh test 2 12 
#
select=$2
ncpus=$3
mpiprocs=$3
ntasks=$((select*mpiprocs))
mpistring="mpirun -np 4"
newmpistring="mpirun -np ${ntasks}"

mkdir -p run
sed -e "s/\${select}/${select}/" -e "s/\${ncpus}/${ncpus}/" -e "s/\${mpiprocs}/${mpiprocs}/" template.pbs > pbs_header.pbs

for f in ./scripts/*.pbs; do 
    nf=$(echo $f | rev | cut -d'/' -f-1 | rev)
    sed -e '1r pbs_header.pbs' -e'/^#/d' $f > ./run/$nf
    sed -i '3i#PBS -N '"$1"'' ./run/$nf 
    sed -i 's/'"$mpistring"'/'"$newmpistring"'/g' ./run/$nf
done

cp ./scripts/PBS_run*.sh ./run/
