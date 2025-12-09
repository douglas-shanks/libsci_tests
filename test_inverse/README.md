A simple benchmark to test performing a matrix inverse using LAPACK.

We first perform an LU decomposition of the matrix A using *getrf, and then 
compute A^{-1} using *getri. 

We then compute x=A^{-1}b using gemv, this allows us to then report the accuracy 
by computing A*x-b = epsilon and printing this value (where epsilon << 0).

Compile: e.g.

ftn -g -homp tester_complex.f90 -o tester_complex

Run: e.g.

time; export OMP_NUM_THREADS=1; srun -N1 -n1 -c${OMP_NUM_THREADS} ./tester_complex

Douglas Shanks 03/12/25
