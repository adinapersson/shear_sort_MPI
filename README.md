# Shear Sort algorithm using MPI

A shear sort algorithm implemented in parallel using MPI in C. The algorithm sorts N = n^2 numbers in a square matrix A of size nxn.

## Run

The program shearsort.c is built by the command "make". 

Run by the command line: "mpirun -np $p shearsort $n input.txt output.txt", where $p is the number of processes and $n is the number of rows/columns in the square matrix. The two last arguments are text files, the input file which the initial data will be read from and the output file which the sorted data will be written to. The input file should first have the number of rows/columns written and then the elements separated by white spaces. 
