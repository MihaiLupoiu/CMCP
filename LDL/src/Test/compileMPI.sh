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

cd ../MPI/
make
mv ./mpi.out ../Test/.

cd ../Test

for (( i=0; i <= $1; i=i+64 ))
do
    for (( j=1; j <= $2; j++ ))
    do
        echo $i
        mpirun -n $j ./mpi.out $i >> tiempoMPI.txt
    done
done

rm *.out
