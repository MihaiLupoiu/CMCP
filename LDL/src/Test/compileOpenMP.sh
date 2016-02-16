#!/usr/bin/env bash

ver=$(uname -s)

if [ -z "$1" ]
then
    echo "First argument not supplied, compile.sh <size>(pow(size>1,2)) <version>(sec,omp,mpi) <cpuNumber>(cpuNumber>500) ";
    exit;
fi

if [ -z "$2" ]
then
    echo "Second argument not supplied, compile.sh <size>(pow(size>1,2)) <cpuNumber>(cpuNumber>500) ";
    exit;
fi

if [ "$ver" == "Darwin" ]; then
    cd ../OpenMP/
    make mac
    mv ./openmp.out ../Test/.

else
    cd ../OpenMP/
    make
    mv ./openmp.out ../Test/.
fi

cd ../Test

for (( i=0; i <= $1; i=i+64 ))
do
    for (( j=1; j <= $2; j++ ))
    do
        echo $i
        ./openmp.out $i $j >> tiempoOpenMP.txt
    done
done

rm *.out
