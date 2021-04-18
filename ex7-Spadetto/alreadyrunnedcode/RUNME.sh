#!/bin/bash                                                                     
echo compiling..                                                  
gfortran ex7.F90 -llapack -lfftw3 -o test     
echo running                                     
./test      
echo gnuplot script                                                                    
gnuplot "gscript.txt"                                                          
echo gnuplotscript2                                                                                
gnuplot "gscript2.txt"                                                                                 
                                                                                
                                                                                
                             
