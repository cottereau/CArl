#!/bin/bash

mkdir -p run
sed -e "s/\${select}/${1}/" -e "s/\${ncpus}/${2}/" -e "s/\${mpiprocs}/${2}/" template.pbs > pbs_header.pbs
for f in ./scripts/*.pbs; do 
    nf=$(echo $f | rev | cut -d'/' -f-1 | rev)
    echo $f 
    echo $nf
    sed -e '1r pbs_header.pbs' -e'/^#/d' $f > ./run/$nf
    sed -i '3i#PBS -N '"$1"'' ./run/$nf 
done

