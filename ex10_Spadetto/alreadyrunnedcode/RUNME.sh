#!/bin/bash
echo "##### STARTING EX10 CODE ######"
rm -r results
mkdir results
gfortran ex10.F90 -llapack -o test 
echo '-code compiled' 
./test

gnuplot g3
