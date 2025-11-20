program complex_matrix_multiply
    implicit none
    
    ! Parameters
    integer :: i,j
    integer, parameter :: N = 1000 ! Size of matrices
    integer, parameter :: lda = N ! Leading dimension of matrices A and C
    integer, parameter :: ldb = N ! Leading dimension of matrix B
    integer, parameter :: ldc = N ! Leading dimension of matrix C
    
    ! Complex numbers
    complex(kind=4) :: alpha, beta
    complex(kind=4), allocatable :: A(:,:), B(:,:), C(:,:)
    real(kind=4) :: rx = 0.0
    real(kind=4) :: ry = 0.0
    real(kind=4) :: output
    real(kind=4), allocatable :: work(:)
    integer :: randomSeedSize

    ! Function declaration for zgemm
    external cgemm
    external clange

    real(kind=4) clange

    allocate(work(N))
    allocate(A(N,N))
    allocate(B(N,N))
    allocate(C(N,N))

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
            A(i,j) = cmplx( rx, ry)
            call random_number(rx)
            call random_number(ry)
            B(i,j) = cmplx( rx, ry)
        end do
    end do
    
    ! Print out matrices A and B
    !print *, "Matrix A:"
    !do i = 1, N
    !    print *, (A(i, j), j = 1, N)
    !end do

    !print *, "Matrix B:"
    !do i = 1, N
    !    print *, (B(i, j), j = 1, N)
    !end do

    print *, ""
    ! Call zgemm to perform matrix multiplication: C = alpha * A * B + beta * C
    alpha = (1.0_4, 0.0_4) ! Real part 1.0, imaginary part 0.0
    beta = (0.0_4, 0.0_4)
    call cgemm('N', 'N', N, N, N, alpha, A, lda, B, ldb, beta, C, ldc)
    
    ! Output the result
    !print *, "Matrix C after multiplication:"
    !do i = 1, N
    !    print *, (C(i, j), j = 1, N)
    !end do

    ! Compute 2 norm of Matrix C
    output = clange('F', N, N, C, N, work)

    print *, "2-norm of Matrix C:"
    print *, output

    deallocate(work, A, B, C)
end program complex_matrix_multiply

