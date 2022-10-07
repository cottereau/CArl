#!/bin/bash
# Run example named "test" on 2 nodes - 24 MPI processes (12 MPI processes per node):
#
#   ./prepare_batch_run.sh pbs test 2 12 
#
if [ $1 == 'pbs' ]
then
    select=$3
    ncpus=$4
    mpiprocs=$4
    ntasks=$((select*mpiprocs))
    mpistring="mpirun -np 4"
    newmpistring="mpirun -np ${ntasks}"

    mkdir -p run
    sed -e "s/\${select}/${select}/" -e "s/\${ncpus}/${ncpus}/" -e "s/\${mpiprocs}/${mpiprocs}/" template.pbs > pbs_header.pbs

    for f in ./scripts/*.pbs; do 
        nf=$(echo $f | rev | cut -d'/' -f-1 | rev)
        sed -e '1r pbs_header.pbs' -e'/^#/d' $f > ./run/$nf
        sed -i '3i#PBS -N '"$2"'' ./run/$nf 
        sed -i 's/'"$mpistring"'/'"$newmpistring"'/g' ./run/$nf
    done

    cp ./scripts/PBS_run*.sh ./run/

elif [ $1 == 'slurm' ] 
then
    select=$3
    ncpus=$4
    mpiprocs=$4
    ntasks=$((select*mpiprocs))
    mpistring="srun -n 4"
    newmpistring="srun -n ${ntasks}"
    ntasksstring="--ntasks=4"
    newntasksstring="--ntasks=${ntasks}"
    nodestring="--nodes=1"
    newnodestring="--nodes=${select}"

    mkdir -p run
    sed -e "s/\${nodes}/${select}/" -e "s/\${ncpus}/${ncpus}/" -e "s/\${mpiprocs}/${mpiprocs}/" template.slurm > slurm_header.slurm

    for f in ./scripts/*.slurm; do 
        nf=$(echo $f | rev | cut -d'/' -f-1 | rev)
        sed -e '1r slurm_header.slurm' -e'/^#/d' $f > ./run/$nf
        sed -i '3i#SBATCH --job-name='"$2"'' ./run/$nf 
        sed -i 's/'"$mpistring"'/'"$newmpistring"'/g' ./run/$nf
        sed -i 's/'"$nodestring"'/'"$newnodestring"'/g' ./run/$nf
        sed -i 's/'"$ntasksstring"'/'"$newntasksstring"'/g' ./run/$nf
    done

    cp ./scripts/SLURM_run*.sh ./run/

fi
