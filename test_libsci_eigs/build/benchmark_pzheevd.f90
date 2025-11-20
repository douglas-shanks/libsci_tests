module helper

implicit none
contains

  elemental subroutine convert_to_int(str, i, stat)
  character(len=*), intent(in) :: str
  integer, intent(out) :: i, stat
  read(str, *, iostat=stat) i
  end subroutine convert_to_int
  
end module helper

program pzheevd_test
  use mpi
  use helper
  implicit none
  integer, parameter :: MY_NPROW = 2
  integer, parameter :: MY_NPCOL = 2
  integer, parameter :: MY_N = 1024
  integer, parameter :: MY_NB = 256
  integer, parameter :: MY_IL = 1
  integer, parameter :: MY_IU = 1024
  integer, parameter :: ELS_TO_PRINT = 5
  
  integer :: info = 0
  character :: jobz = 'V'
  character :: uplo = 'L'
  integer :: ln = MY_N
  integer :: lnbrow = MY_NB
  integer :: lnbcol
  integer :: m, nz
  integer :: err
  integer :: il = MY_IL, iu = MY_IU
  real(kind=8) :: dummyL, dummyU
  real(kind=8) :: t1, t2
  integer :: lnprow = MY_NPROW
  integer :: lnpcol = MY_NPCOL
  integer :: myrow, mycol
  integer :: ia=1, ja=1, iz=1, jz=1
  integer :: desca(15), descz(15)
  integer :: ctxt_sys
  integer :: moneI = -1, zeroI = 0, oneI = 1
  integer :: rank, size, i
  integer :: llda
  character(len=9) :: procOrder = "Row-major"
  character(len=10), dimension(4) :: argc

  ! ADDED DEFINITIONS
  integer :: ml, nl
  integer :: iseed(4)
  integer :: numEle
  complex(kind=8), allocatable :: a_ref(:,:), z(:,:)
  real(kind=8), allocatable :: w(:)
  complex(kind=8), allocatable :: work(:)
  real(kind=8), allocatable :: rwork(:)
  integer, allocatable :: iwork(:)
  real(kind=8) :: temp(2), rtemp(2)
  integer :: liwork, lwork, lrwork
  integer :: my_rank, j, ierror
  
  ! .. External Functions ..
  integer, external :: numroc

  ! Initialize MPI and BLACS
  call MPI_INIT(err)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, err)
  do i = 1, 4
    call get_command_argument(i, argc(i))
  end do

  call convert_to_int(argc(1), ln, err)
  call convert_to_int(argc(2), lnprow, err)
  call convert_to_int(argc(3), lnpcol, err)
  call convert_to_int(argc(4), lnbrow, err)
  lnbcol = lnbrow

  if (rank == 0) then
    print *, "ln =", ln
    print *, "lnprow =", lnprow
    print *, "lnpcol =", lnpcol
    print *, "lnbropw =", lnbrow
  end if
  
  call blacs_get(moneI, zeroI, ctxt_sys)
  call blacs_gridinit(ctxt_sys, procOrder, lnprow, lnpcol)
  call blacs_gridinfo(ctxt_sys, lnprow, lnpcol, myrow, mycol)
  
  if (myrow .eq. -1) then
    print *, "Failed to initialize MPI and/or BLACS!"
    call MPI_FINALIZE(err)
    stop
  end if

  ! Allocate my matrices now
  ml = numroc(ln, lnbrow, myrow, zeroI, lnprow)
  nl = numroc(ln, lnbrow, mycol, zeroI, lnpcol)
  llda = ml
  call descinit(desca, ln, ln, lnbrow, lnbrow, &
                & zeroI, zeroI, ctxt_sys, llda, info)
  call descinit(descz, ln, ln, lnbrow, lnbrow, &
                & zeroI, zeroI, ctxt_sys, llda, info)
  allocate( a_ref(ml,nl) )
  allocate( z(ml,nl) )
  allocate( w(ln) )
  iseed(1) = myrow
  iseed(2) = mycol
  iseed(3) = mycol + myrow*lnpcol
  iseed(4) = 1
  
  if (iand(iseed(4), 2) == 0) then
    iseed(4) = iseed(4) + 1
  end if
  
  numEle = ml*nl
  call zlarnv(oneI, iseed, numEle, a_ref)
  call zlarnv(oneI, iseed, numEle, z)
  t1 = MPI_WTIME()
  call pzheevd(jobz, uplo, ln, a_ref, ia, ja, &
        & desca, w, z, iz, jz, descz, temp, moneI, &
        & rtemp, moneI, liwork, moneI, info)
  lwork = temp(1)
  
  allocate( work(lwork) )
  lrwork = rtemp(1)
  allocate( rwork(lrwork) )
  allocate( iwork(liwork) )
  
  call pzheevd(jobz, uplo, ln, a_ref, ia, &
        & ja, desca, w, z, iz, jz, descz, work, &
        & lwork, rwork, lrwork, iwork, liwork, info)
        
  if (info /= 0) then
    write(6,*) "PZHEEVD returned non-zero info val of ", info
  end if
  
  t2 = MPI_WTIME()
  
  if (rank == 0) then
    print *, "PZHEEVD time: ", t2 - t1
  end if
  
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierror)
  call blacs_gridexit(ctxt_sys)
  call blacs_exit(zeroI)
  
  ! if (rank == 0 ) then
  !   do j=1, ELS_TO_PRINT
  !   write(6,*) "w(", j, ")=", w(j), "z(", j, ")=", z(j,1)
  !   end do
  ! end if
  
  deallocate( a_ref, z, w )
  deallocate( work, rwork, iwork )
  
end program pzheevd_test