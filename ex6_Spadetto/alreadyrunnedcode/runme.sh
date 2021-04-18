#!/bin/bash 
gfortran harmonic.F90 -llapack -o harm -fcheck=all -Wall
./harm
cd results
gnuplot "../gnuplotscripts/gscript"
gnuplot "../gnuplotscripts/gscript2"
gnuplot "../gnuplotscripts/gscript3"
gnuplot "../gnuplotscripts/gscript5"
