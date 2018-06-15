!
   program test_mpdgeqrf

!-------------------------------------------------------------------------------
! General matrix QR factorization - DOUBLE PRECISION version
!-------------------------------------------------------------------------------

   implicit none
   include 'mpif.h'

   !-------------------------------------------------------------------------------
   ! Please set system size parameters below!
   ! ma:   matrix row dimension
   ! na:   matrix column dimension 
   ! nblk: Blocking factor in block cyclic distribution
   !-------------------------------------------------------------------------------

   integer :: mblk, nblk
   integer :: ma, na

   !-------------------------------------------------------------------------------
   !  Local Variables

   integer np_rows, np_cols, na_rows, na_cols, lwork, &
           IID, JJD, ID, IDROW, IDCOL

   integer myid, nprocs, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols
   integer i, mpierr, my_blacs_ctxt, sc_desc(9), info, nprow, npcol,itmp(4)

   integer, external :: numroc
   integer, allocatable :: iwork(:)

   real*8 err, errmax, ttt0, ttt1
   real*8, allocatable :: ev(:)
   real*8 :: xc,rtmp(4)
   real*8, allocatable :: a(:,:), as(:,:), work(:), tmp1(:), tmp2(:)

   real*8, parameter :: ZERO = 0.d0, ONE = 1.d0, NEONE=-1.0D+00

   integer :: iseed(4096)     ! Random seed, size should be sufficient for every generator

   integer :: STATUS
   real*8, external  :: pdlange
   !-------------------------------------------------------------------------------
   !  Parse command line argumnents, if given

   mblk = 64 
   nblk = 32 
   ma  = 10000
   na = 5000
 
   !-------------------------------------------------------------------------------
   ! MPI Initialization
   call mpi_init( mpierr )
   call mpi_comm_rank( mpi_comm_world,myid,mpierr )
   call mpi_comm_size( mpi_comm_world,nprocs,mpierr )

   STATUS = 0

!  Initialize the BLACS: refer to user's guide or lbll manual
   call BLACS_PINFO( myid,nprocs )

   !-------------------------------------------------------------------------------
   ! Selection of number of processor rows/columns
   ! We try to set up the grid square-like, i.e. start the search for possible
   ! divisors of nprocs with a number next to the square root of nprocs
   ! and decrement it until a divisor is found.

   do np_rows = NINT(SQRT(REAL(nprocs))),2,-1
      if(mod(nprocs,np_rows) == 0 ) exit
   enddo
   ! at the end of the above loop, nprocs is always divisible by np_cols

   np_cols = nprocs/np_rows

   if(myid==0) then
      print *
      print '(a)','Standard eigenvalue problem - DOUBLE PRECISION version'
      print *
      print '(3(a,i0))','Matrix size=',ma,' -by- ',na,', Block size=',nblk
      print '(3(a,i0))','Number of processor rows=',np_rows,', cols=',np_cols,', total=',nprocs
      print '(a,i0)', 'mpi_comm_world, ', mpi_comm_world
      print *
   endif

   !-------------------------------------------------------------------------------
   ! Set up BLACS context and MPI communicators
   !
   ! The BLACS context is only necessary for using Scalapack.
   !
   ! For ELPA, the MPI communicators along rows/cols are sufficient,
   ! and the grid setup may be done in an arbitrary way as long as it is
   ! consistent (i.e. 0<=my_prow<np_rows, 0<=my_pcol<np_cols and every
   ! process has a unique (my_prow,my_pcol) pair).

!  Initialize a single BLACS context

   my_blacs_ctxt = mpi_comm_world
   call BLACS_GET( -1,0,my_blacs_ctxt )
   call BLACS_Gridinit( my_blacs_ctxt, 'c', np_rows, np_cols )
   call BLACS_Gridinfo( my_blacs_ctxt, nprow, npcol, my_prow, my_pcol )

   ! Determine the necessary size of the distributed matrices,
   ! we use the Scalapack tools routine NUMROC for that.

   na_rows = numroc(ma, mblk, my_prow, 0, np_rows)
   na_cols = numroc(na, nblk, my_pcol, 0, np_cols)
   ! write(*,*) 'myprow, mypcol', my_prow, my_pcol, na_rows, na_cols

   ! Set up a scalapack descriptor for the checks below.
   ! For ELPA the following restrictions hold:
   ! - block sizes in both directions must be identical (args 4+5)
   ! - first row and column of the distributed matrix must be on row/col 0/0 (args 6+7)

   call descinit( sc_desc, ma, na, mblk, nblk, 0, 0, my_blacs_ctxt, na_rows, info )
   if(myid == 0 ) then 
     print *
     write(*,*) 'sc_desc', sc_desc(1:9)
     write(*,*) 'narow and nacols', na_rows, na_cols
     print *
   endif

   !-------------------------------------------------------------------------------
   ! Allocate matrices and set up a test matrix for the eigenvalue problem
   allocate( a (na_rows,na_cols) )
   allocate( as(na_rows,na_cols) )
   allocate( ev(na) )

   iseed(:) = myid
   call RANDOM_SEED(put=iseed)
   call RANDOM_NUMBER(a)

   CALL PDLASET( 'FULL',na,na,ZERO,ZERO,a,na+1,1,sc_desc )
   ! Append the identity matrix
   DO ID=1, na, 1
      CALL INFOG2L( ID+na,ID,sc_desc,np_rows,np_cols,my_prow,my_pcol,&
           IID,JJD,IDROW,IDCOL )
      IF( my_prow.EQ.IDROW .AND. my_pcol.EQ.IDCOL ) THEN
         A( IID,JJD ) = ONE
      END IF
   END DO

   ! Save original matrix A for later accuracy checks
   as = a
   !-------------------------------------------------------------------------------

   !-------------------------------------------------------------------------------
   ! Calculate eigenvalues/eigenvectors by calling pzheevd
   lwork=-1
   rtmp(:) =ZERO
   itmp(:) =0

   call MPDGEQRF( Ma, Na, A, 1, 1, sc_desc, ev, rtmp, LWORK, INFO )

   if( info .ne. 0 ) then
      print *, 'pdgeqrf returns error during 1st call', info
   end if

   lwork=2*int(rtmp(1))+1

   IF( myid==0 ) print*, 'lwork, ', lwork

   allocate( work(lwork),tmp1(na_rows*na_cols), tmp2( na_rows*na_cols ) )

   ttt0 = MPI_Wtime()
   call MPDGEQRF( Ma, Na, A, 1, 1, sc_desc, ev, work, LWORK, INFO )
   ttt1 = MPI_Wtime()
   if(myid == 0) print *,'QR factorization costs  :',ttt1-ttt0, 'info', info

   ! Copy R to tmp1
   CALL PDLASET( 'FULL',ma,na,ZERO,ZERO,tmp1,1,1,sc_desc )
   CALL PDLACPY( 'Up',na,na,A,1,1,sc_desc,tmp1,1,1,sc_desc ) 

   ttt0 = MPI_Wtime()
   call MPDORGQR( Ma, Na, Na, a, 1, 1, sc_desc, ev, WORK, LWORK, INFO )
   ttt1 = MPI_Wtime()
   IF( myid == 0) print *,'Generate Q explicitly costs:',ttt1-ttt0, 'info', info

   CALL PDGEMM( 'N','N',Ma,Na,Na,ONE,A,1,1,sc_desc,tmp1,1,1,sc_desc,ZERO,&
                 tmp2,1,1,sc_desc  ) 
 
   ! CALL PDLAPRNT( ma,na,tmp2,1,1,sc_desc,1,1,'B1',6,WORK )  

   call pdgeadd( 'N',Ma,Na,NEONE,As,1,1,sc_desc,ONE,tmp2,1,1,sc_desc )
   err = PDLANGE( 'Fro',Ma,Na,tmp2,1,1,sc_desc,WORK )

   IF( myid ==0 ) write(*,*) 'err= ', err
  
   deallocate(a)
   deallocate(as)
   deallocate(ev)
   deallocate(tmp1)
   deallocate(tmp2)
   deallocate(work)

   call blacs_gridexit(my_blacs_ctxt)
!   call mpi_finalize(mpierr)
   call EXIT(STATUS)

   end
