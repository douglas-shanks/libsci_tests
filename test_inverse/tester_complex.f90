Program test_inv

    implicit none

    integer, parameter :: dp = kind(1.d0)

    real(kind=dp) :: dznrm2
    external dznrm2
    
    complex(kind=dp), dimension(:,:), allocatable :: A,B
    complex(kind=dp), dimension(:), allocatable :: WORK,C,D,Y,X
    real(kind=dp) :: norm,div
    integer, dimension(:), allocatable :: IPIV
    integer :: i, j, SIZ, INFO

    SIZ=4000

    allocate(A(SIZ,SIZ))
    allocate(B(SIZ,SIZ))
    allocate(C(SIZ))
    allocate(D(SIZ))
    allocate(X(SIZ))
    allocate(Y(SIZ))
    allocate(WORK(SIZ))
    allocate(IPIV(SIZ))

    div=real(1.0/2.0,8)
    do i=1, SIZ
        do j=1, SIZ
            A(i, j)=cmplx(ranf(),ranf())*div+div
            B(i, j)=0.0_8
            C(i)=cmplx(ranf(),ranf())*div+div
            D(i)=0.0_8
        end do
        !ensure diagonal dominance
        A(i,i)=A(i,i)+1.0_8
    end do

    B=A
    Y=D

    !print *, A
    call zgetrf(size(B, 1), size(B, 2), B, size(B, 1), IPIV, INFO)
    !print *, "========="
    !print *, B
    !print *, IPIV
    call zgetri(size(B, 1), B, size(B, 1), IPIV, WORK, size(WORK), INFO);
    !print *, "========="
    !print *, B
    

    !Solve A

    !
    ! D=A^{-1}C
    ! DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
    call zgemv('N',size(B, 1), size(B, 2),1.0_8,B,size(B, 1),C,1,0.0_8,D,1)
    !print*, "D", D
    ! A*D-C=0
    call zgemv('N',size(B, 1), size(B, 2),1.0_8,A,size(B, 1),D,1,0.0_8,Y,1) 
    !print*, "Y", Y
    X=Y-C
    !print*, "X", X
    norm = dznrm2(size(X,1),X,1)
    print*, "size(A)", size(A)
    print*, "L2 norm of Ax-b", norm
    

    deallocate(A)
    deallocate(B)
    deallocate(C)
    deallocate(D)
    deallocate(X)
    deallocate(Y)
    deallocate(WORK)
    deallocate(IPIV)

end program test_inv
