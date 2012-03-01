#!/bin/bash
mkdir $1
gcc -O1 -g -lm main/prepare_job.c -o $1/prepare_job
cd $1
cp ../lasagna* .
cp ../*.dat .
./prepare_job $1
chmod +x runlasagna*
chmod +x master.sh
