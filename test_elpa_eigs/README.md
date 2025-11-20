A simple benchmark to test ELPA eigenvalue solvers.

Adapted from code in Appendix of paper "Eigensolver Performance Comparison on Cray XC Systems" by B. Cook et al at CUG 2018.

Currently only PDSYEVD takes input from input.dat which specifies for a NxN real matrix
nrow, ncol, N, nb

Where:
  -  nrow  = MPI ranks in x
  -  nrcol = MPI ranks in y
  -  nb    = block size.

Remember that:
  - total MPI ranks, nmpi=nrow x ncol and
  - N = nrow x ncol x nb

Douglas Shanks 07/11/23
