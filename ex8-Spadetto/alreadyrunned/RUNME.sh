#!/bin/bash                                                                     
echo compiling..                                                  
gfortran ex8.F90 -llapack -lfftw3 -o test     
echo running                                     
./test      
cd results
echo gnuplot script                                                                    
gnuplot "../text"                                                          
                                                                  
                                                                                
                                                                                
                             
