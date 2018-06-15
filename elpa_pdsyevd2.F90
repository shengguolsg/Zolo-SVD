!
  SUBROUTINE ELPA_PDSYEVD2( UPLO,na,nev,A,LDA,LCA,EV,Z,nblk,sc_desc,np_rows,&
                np_cols,my_prow,my_pcol,myid,mpi_comm_In,info,REAL_KERNEL ) 
!
!
    use ELPA1
    use ELPA2 
    use test_util

    implicit none
    include 'mpif.h'
!
! .. Scalar Arguments ..
    CHARACTER, INTENT(IN)  :: UPLO
    INTEGER, INTENT(IN)    :: na,nev,LDA,LCA,nblk,my_prow,my_pcol,myid, & 
                              np_rows,np_cols
    INTEGER, INTENT(INOUT) :: info
    INTEGER, optional  :: mpi_comm_In, REAL_KERNEL
!
! .. Array Arguments ..
    REAL*8, INTENT(INOUT)     :: A(LDA,*), Z(LDA,*)
    REAL*8, INTENT(INOUT)     :: EV(*)
    INTEGER, INTENT(IN)       :: sc_desc(*)
!
! .. Purpose ..
! ==============
!
!  Solves a real symmetric eigenvalue problem with a 2 stage approach
!                A* Q = \lambda Q,                      (*)
!  where A is hermitian. The purpose of this routine is to 
!  write an equivalent interface to pzheevx, and therefore we can easily
!  replace
!  pzheevx with elpa routine. 
!
!  Note that we need to initialize MPI communicator before calling elpa
!  routines since elpa needs the row and column MPI communicators. 
!
! ..Parameters..
! ===============
!
!  UPLO    (global input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          real symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  na          Order of matrix A
!
!  nev         Number of eigenvalues needed
!
!  A(LDA,*)    Distributed matrix for which eigenvalues are to be
!  computed.
!              Distribution is like in Scalapack.
!              The full matrix must be set (not only one half like in
!              scalapack).
!              Destroyed on exit (upper and lower half).
!
!  LDA         Leading dimension of A, local dimension. A, B and Z have
!  the same
!              leading dimension. 
!
!  EV(na)      On output: eigenvalues of a, every processor gets the
!  complete set
!
!  Z(LDA,*)    On output: Eigenvectors of a
!              Distribution is like in Scalapack.
!              Must be always dimensioned to the full size
!              (corresponding to (na,na))
!              even if only a part of the eigenvalues is needed.
!
!  nblk        blocksize for 2D cyclic distribution, must be the same in
!  both directions!
!
!  sc_desc     An scalapack descriptor, an integer array with
!  dimension(9), desribe 
!              the data distribution.
!
!  np_rows     the row dimension of current process grid 
!
!  np_cols     the col dimension of current process grid 
!
!  my_prow     the row index of current process in the process grid
!
!  my_pcol     the column index of current process in the process grid
!
!  myid        the ID of current processor
! 
!  info        integer output, signal success if equal 0 or failure
!
!  mpi_comm_In 
!              Optional parameter, If present, it equals the mpi
!              Communicator of current processors
!                                  If not,     it equals all
!                                  MPI_COMM_WORLD by default
!              It is usefull for multi level parallelization case. 
!  REAL_KERNEL
!              Optional parameter, choosing a real kernel for 2-stage
!              method
!
! =================================================================
! Written by Shengguo Li, in Aug. 19, 2014
! Another interface of solve_evp_real_2stage for the standard eigenvalue
! problem.
!
! If you meet problem and are familiar with ELPA, you can try to use
! solve_evp_real_2stage.
! =================================================================
!
! .. Local Scalars ..
! 
    integer na_rows, i, n_col, n_row, na_cols, use_real2_kernel
    integer mpi_comm_rows,mpi_comm_cols,mpi_comm_all, mpierr
    logical success, lower, upper, full
    real*8  ttt0, ttt1
    real*8, parameter :: ZERO = 0.d0, ONE = 1.d0

    real*8, allocatable :: WORK(:,:)

    integer, external :: numroc, indxl2g
    logical, external :: lsame
!
    info = 0
    success = .true. 
    na_rows = LDA
    na_cols = LCA

! ****** Check for input parameters *******
   lower = lsame(UPLO,'L')
   upper = lsame(UPLO,'U')
   full  = lsame(UPLO,'A')

   IF (.NOT. (lower .OR. upper .OR. full) ) Then
      print *, 'The first parameter of elpa_pdsyevd must be L or U or A'
      info = -1
      return
   ENDIF

!  ELPA solver needs a full matrix. If the input matrix A is triangular
!  form,
!  set the other half part. In practice, A is usually constructed fully. 
   IF (.NOT. full ) THEN 
  
      allocate( work(LDA,LCA) )
      call pdtran(na,na,ONE,a,1,1,sc_desc,ZERO,work,1,1,sc_desc)
 
      IF (upper) THEN
        do i=1,na_cols
           ! Get global column corresponding to i and number of local
           ! rows up to
           ! and including the diagonal, these are unchanged in A
           n_col = indxl2g(i,     nblk, my_pcol, 0, np_cols)
           n_row = numroc (n_col, nblk, my_prow, 0, np_rows)
           a(n_row+1:na_rows,i) = work(n_row+1:na_rows,i)
       enddo
      ELSE IF (lower ) THEN
        do i=1,na_cols
           ! Get global column corresponding to i and number of local
           ! rows up to
           ! and including the diagonal, these are unchanged in A
           n_col = indxl2g(i,     nblk, my_pcol, 0, np_cols)
           n_row = numroc (n_col, nblk, my_prow, 0, np_rows)
           a(1:n_row-1,i) = work(1:n_row-1,i)
       enddo
     END IF

     deallocate( work )
   END IF

    ! set print flag in elpa1
    elpa_print_times = .false.

#ifdef Output
    elpa_print_times = .true.
    IF(myid ==0 ) print *, 'Start to use elpa Two stage...'
#endif

   ! All ELPA routines need MPI communicators for communicating within
   ! rows or columns of processes, these are set in
   ! get_elpa_row_col_comms.
   IF ( present(mpi_comm_In) ) THEN
       mpi_comm_all = mpi_comm_In
   ELSE
       mpi_comm_all = mpi_comm_world
   END IF
   mpierr= get_elpa_communicators( mpi_comm_all, my_prow, my_pcol, &
                               mpi_comm_rows, mpi_comm_cols )

   ! use_real2_kernel can equal (1,2),5,6,7,8. It is better not to
   ! choose 
   ! (1,2). You can try the last four choices.    
   IF ( present(REAL_KERNEL) ) THEN
       use_real2_kernel = REAL_KERNEL 
   ELSE
       use_real2_kernel = 12 
   END IF

   call mpi_barrier(mpi_comm_all, mpierr) ! for correct timings only
   ttt0 = MPI_Wtime()
   success = solve_evp_real_2stage(na, nev, a, na_rows, ev, z, na_rows,nblk, &
                 na_cols,mpi_comm_rows, mpi_comm_cols, mpi_comm_all, use_real2_kernel )
   ttt1 = MPI_Wtime()
   if (myid == 0 ) write(*,*) 'elpa2 costs ', ttt1 - ttt0
  
   if (.not.(success)) then
      write(*,*) "solve_evp_real_2stage in elpa_pdsyevd2 produced an error! Aborting..."
      info = 9 
      call MPI_ABORT(mpi_comm_all, mpierr)
   endif

#ifdef Output
   if(myid == 0) print *,'Time transform to tridi :',time_evp_fwd
   if(myid == 0) print *,'Time solve tridi        :',time_evp_solve
   if(myid == 0) print *,'Time transform back EVs :',time_evp_back
   if(myid == 0) print *,'Total time (sum above)
:',time_evp_back+time_evp_solve+time_evp_fwd
#endif


 end subroutine elpa_pdsyevd2
