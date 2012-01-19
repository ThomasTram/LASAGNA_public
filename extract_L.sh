#!/bin/bash
mkdir output
plotstr=plot
for i in {0..2}
do
  for j in {0..2}
  do
	filename=dump_00${i}_00${j}.mat
	../extract_matrix ${filename} T L
	paste output/T.dat output/L.dat > output/L_${i}_${j}.dat
	plotstr="$plotstr 'output/L_${i}_${j}.dat' u 1:2 w l",
  done
done
plotstr="$plotstr 0 0"
echo "set logscale y" >script.gnu
echo ${plotstr} >> script.gnu
