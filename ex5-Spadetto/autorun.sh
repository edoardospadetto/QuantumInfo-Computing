#!/bin/bash
gcc -E test5.F90 -o _test5.F90
gfortran _test5.F90 -llapack -o ./test
./test
cd results
gnuplot "g_script.txt"
