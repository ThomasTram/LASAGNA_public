#!/bin/bash
mkdir $1
cp main/prepare_job.c $1
cd $1
gcc -O1 -g -lm prepare_job.c -o prepare_job
cp ../lasagna* .
cp ../*.dat .
./prepare_job $1
chmod +x runlasagna*
chmod +x master.sh
