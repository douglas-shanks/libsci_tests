program real_matrix_solve
    implicit none
    
    ! Parameters
    integer :: i,j
    integer, parameter :: N = 1000 ! Size of matrices
    integer, parameter :: lda = N ! Leading dimension of matrices A and C
    integer, parameter :: ldb = N ! Leading dimension of matrix B
    
    ! Real numbers for solving Ax = b
    real(kind=8), allocatable :: A(:,:), As(:,:)
    real(kind=8), allocatable :: b(:), x(:), y(:), diff(:)
    real(kind=8), allocatable :: pivot(:) ! Pivot indices
    real(kind=8) :: rx = 0.0
    real(kind=8) :: ry = 0.0
    real(kind=8) :: norm
    real(kind=8), allocatable :: work(:)
    integer :: randomSeedSize
    integer :: info

    ! Function declaration for dgesv
    external dgesv, dgemv, scnrm2

    real(kind=8) scnrm2

    allocate(A(N,N), As(N,N), b(N), x(N), y(N), diff(N), pivot(N))

    ! Initialize matrices A and B
    ! Seed random number generator
    call random_seed(size = randomSeedSize)
    call random_seed()  ! use system-provided seed

    ! Generate random complex matrices A and B
    do i = 1, N
        do j = 1, N
            ! Generate random real and imaginary parts between 0 and 1
            call random_number(rx)
            call random_number(ry)
            A(i,j) = rx
            b(i) = ry
        end do
    end do
            
    As=A
    x=b
    y=0.0_8

    ! Print out matrices A and B
    !print *, "Matrix A:"
    !do i = 1, N
    !    print *, (A(i, j), j = 1, N)
    !end do

    !print *, "Vector B:"
    !do i = 1, N
    !    print *, b(i)
    !end do

    ! Call dgesv to perform matrix solve. Solution overwritten on b
    call dgesv(N, 1, A, lda, pivot, b, ldb, info)
    
    ! Output the result
    !print *, "Solution:"
    !do i = 1, N
    !    print *, b(i)
    !end do

    ! Check solution A*x = b

    call dgemv('N', N, N, 1.0_8, As, lda, b, 1, 0.0_8, y, 1)

    ! Output the result
    !print *, "Solution:"
    !do i = 1, N
    !    print *, y(i)
    !end do

    diff = x-y
    norm = scnrm2(N,diff,1)

    if (norm <= 1.0e-15) then
            print*, "SUCCESS"
    endif 
    !print *, "Norm of difference (should be zero):"
    !print*, norm



    deallocate(A, b, x, y, pivot)

end program real_matrix_solve

