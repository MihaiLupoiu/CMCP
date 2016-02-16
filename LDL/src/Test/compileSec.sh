#!/usr/bin/env bash

ver=$(uname -s)

if [ -z "$1" ]
then
    echo "First argument not supplied, compile.sh <size>(pow(size>1,2))";
    exit;
fi

cd ../Secuencial/
make
mv ./sequential.out ../Test/.

cd ../Test

for (( i=0; i <= $1; i=i+64 ))
do
    echo $i
    ./sequential.out $i >> tiempoSecuencial.txt

done

rm *.out
