The matrixqic.f90 file contains all the function requested.
While programming i wrote more functions than necessary (to play with lapack :) ), so just correct the ones introduced in the .pdf.

subroutine Init
subroutine cgadj 
subroutine ctrace
subroutine ofile
function trace
function adjoint 
the two interfaces
and the type declaration

to compile: 

gfortran -g test_ex21.f90  -o test  -llapack -fcheck=all

to run: 

./test
