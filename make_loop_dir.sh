#!/bin/bash
mkdir $1
gcc -O1 -lm main/prepare_job.c -o $1/prepare_job
cd $1
cp ../lasagna .
cp ../*.dat .
./prepare_job
chmod +x runlasagna*
chmod +x master.sh
