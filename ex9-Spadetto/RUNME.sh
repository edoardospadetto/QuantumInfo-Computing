#!/bin/bash
echo "##### STARTING EX9 CODE ######"
echo "removing old results" 
rm -r results3
rm -r results4
rm -r results5
rm -r results8
rm -r lowest_eig

mkdir results3
mkdir results4
mkdir results5
mkdir results8
echo '-generated directories to store results'
gfortran ex9.F90 -llapack -o test 
echo '-code compiled' 
./test
echo 'code runned'
mkdir lowest_eig
cp results3/data1.txt lowest_eig/data3.txt 
cp results4/data1.txt lowest_eig/data4.txt 
cp results5/data1.txt lowest_eig/data5txt 
cp results8/data1.txt lowest_eig/data8.txt 
 

gnuplot gnuplotscripts/gscript
gnuplot gnuplotscripts/tt
gnuplot gnuplotscripts/g3
