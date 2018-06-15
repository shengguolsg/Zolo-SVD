!
   program pqdwheig

!-------------------------------------------------------------------------------
! This routine first calls pqdwhfac.f90 to compute a polar factorization of an symmetric
! n x n real matrix A. It uses QDWHEIG algorithm to reduce the order of interested eigenvalue
! part. 
!      
!-------------------------------------------------------------------------------
   use mod_prepare_matrix1

   implicit none
   include 'mpif.h'

   !-------------------------------------------------------------------------------
   ! Please set system size parameters below!
   ! na:  the size of matrix A
   ! nq:  the size of matrix Q, nq = 2*na  
   ! nblk: Blocking factor in block cyclic distribution
   !-------------------------------------------------------------------------------

   integer :: nblk
   integer na, nq
!     .. Parameters ..
      INTEGER            BLOCK_CYCLIC_2D, CSRC_, CTXT_, DLEN_, DTYPE_,  &
                         LLD_, MB_, M_, NB_, N_, RSRC_
      PARAMETER          ( BLOCK_CYCLIC_2D = 1, DLEN_ = 9, DTYPE_ = 1,  &
                           CTXT_ = 2, M_ = 3, N_ = 4, MB_ = 5, NB_ = 6, &
                           RSRC_ = 7, CSRC_ = 8, LLD_ = 9 )

   !-------------------------------------------------------------------------------
   !  Local Variables

   integer np_rows, np_cols, na_rows, na_cols, nq_rows, nq_cols, nbw, nnz

   integer myid, nprocs, my_prow, my_pcol, mpi_comm_rows, mpi_comm_cols
   integer i, mpierr, my_blacs_ctxt, desca(9), descw(9),info, nprow, npcol, &
           lwork,iptau,ipwork2, lwork2, lworkt, Rk,liwork, lwork1

   integer, external :: numroc
   character(len=128)   :: filename

   real*8  err, errmax, ttt0, ttt1,norma, sigma
   real*8  rrtmp(4)
   real*8  rtmp(4)
   integer itmp(4)
   real*8, allocatable :: a(:,:), q(:,:), as(:,:), work(:)
   integer,allocatable :: iwork(:)

   real*8, parameter :: ZERO = 0.d0, ONE = 1.d0, NEONE=-1.0D+0

   integer :: iseed(4096)     ! Random seed, size should be sufficient for every generator
   INTEGER, PARAMETER :: nnzmax= 147631 
   INTEGER :: indx(nnzmax), jndx(nnzmax)
   REAL*8  :: rval(nnzmax)

   DOUBLE PRECISION, EXTERNAL :: PDLANGE

   integer :: STATUS
   !-------------------------------------------------------------------------------
   !  Parse command line argumnents, if given

   nblk = 32 
   na  = 153
   nq  = 2*na
   sigma = 1.0D+07
   !-------------------------------------------------------------------------------
   !  MPI Initialization
!   call mpi_init(mpierr)
!   call mpi_comm_rank(mpi_comm_world,myid,mpierr)
!   call mpi_comm_size(mpi_comm_world,nprocs,mpierr)

   STATUS = 0

!  Initialize the BLACS: refer to user's guide or lbll manual
   call BLACS_PINFO(myid,nprocs)

   !-------------------------------------------------------------------------------
   ! Selection of number of processor rows/columns
   ! We try to set up the grid square-like, i.e. start the search for possible
   ! divisors of nprocs with a number next to the square root of nprocs
   ! and decrement it until a divisor is found.

   do np_cols = NINT(SQRT(REAL(nprocs))),2,-1
      if(mod(nprocs,np_cols) == 0 ) exit
   enddo
   ! at the end of the above loop, nprocs is always divisible by np_cols

   np_rows = nprocs/np_cols

   if( myid==0 ) then
      print *
      print '(a)','Test polar factorization - DOUBLE PRECISION version'
      print *
      print '(2(a,i0))','Matrix size=',na, 'Block size=',nblk
      print '(3(a,i0))','Number of processor rows=',np_rows,', cols=',np_cols,', total=',nprocs
      print *
   endif

   !-------------------------------------------------------------------------------
   ! Set up BLACS context and MPI communicators
   !
   ! Initialize a single BLACS context

   !my_blacs_ctxt = mpi_comm_world
   call BLACS_GET(-1,0,my_blacs_ctxt)
   call BLACS_Gridinit( my_blacs_ctxt, 'R', np_rows, np_cols )
   call BLACS_Gridinfo( my_blacs_ctxt, nprow, npcol, my_prow, my_pcol )

   if(myid == 0 ) then
      print '(2(a,i3))', 'np_rows=', np_rows, 'nprow=', nprow
      print '(2(a,i3))', 'np_cols=', np_cols, 'nprow=', npcol
   end if

   ! Determine the necessary size of the distributed matrices,
   ! we use the Scalapack tools routine NUMROC for that.

   na_rows = numroc(na, nblk, my_prow, 0, np_rows)
   na_cols = numroc(na, nblk, my_pcol, 0, np_cols)
   nq_rows = numroc(nq, nblk, my_prow, 0, np_rows)
   nq_cols = numroc(na, nblk, my_pcol, 0, np_cols)

!   write(*,*) 'narows', na_rows, na_cols, my_prow, my_pcol, myid
!   write(*,*) 'nqrows', nq_rows, nq_cols, my_prow, my_pcol, myid

   ! Set up a scalapack descriptor for the checks below.
   ! For ELPA the following restrictions hold:
   ! - block sizes in both directions must be identical (args 4+5)
   ! - first row and column of the distributed matrix must be on row/col 0/0 (args 6+7)

   call descinit( desca, na, na, nblk, nblk, 0, 0, my_blacs_ctxt, na_rows, info )
   call descinit( descw, nq, na, nblk, nblk, 0, 0, my_blacs_ctxt, nq_rows, info )

   !-------------------------------------------------------------------------------
   ! Allocate matrices and set up a test matrix for the eigenvalue problem

   allocate( a (na_rows,na_cols) )
   allocate( q (na_rows,na_cols) )
   allocate( as(na_rows,na_cols) )

   ! For getting a hermitian test matrix A we get a random matrix Z
   ! and calculate A = Z + Z**H

   ! We want different random numbers on every process
   ! (otherways A might get rank deficient):

!   iseed(:) = myid
!   call RANDOM_SEED(put=iseed)
!   call RANDOM_NUMBER(q)

!   a(:,:) = q(:,:)
!   call pdtran(na, na, ONE, q, 1, 1, desca, ONE, a, 1, 1, desca) ! A = Z + Z**T
   
   !filename = 'linverse.mtx'
   filename = 'bcsstk05.mtx'
   CALL PDLASET( 'A',na,na,ZERO,ZERO,A,1,1,desca )
   CALL PDLASET( 'A',na,na,ZERO,ZERO,q,1,1,desca )
   CALL prepare_real_symm_MM( filename,na,nbw,nnz,myid,my_prow,my_pcol, &
        np_rows,np_cols,desca,A,Q,zero )

   ! Save original matrix A for later accuracy checks
   as = a

   !-------------------------------------------------------------------------------
   ! set up workspace and use pqdwh to compute the polar factor of A 
   ! The problem left is to estimate the size of WORK
   LWORK = -1
   call pqdwhfacs( na,na,zero,a,1,1,desca,q,1,1,desca,rtmp,lwork,descw,rrtmp,rtmp,lwork2,info,'N' )
   LWORK = int(rtmp(1)) + 1 + 8*na
   LWORK2 = int(rtmp(2))+1
   LWORKT= LWORK + LWORK2 + na

!   WRITE(*,*) 'The total workspace', lworkt, lwork, lwork2
   
!   LWORK = nq_rows*nq_cols+4*na + na_rows*na_cols

!   write(*,*) 'nqrows', nq_rows, nq_cols, lwork, my_prow, my_pcol
!   write(*,*) 'narows', na_rows, na_cols, my_prow, my_pcol

   allocate( work(lworkt) )
   IPTAU = lwork +1 
   IPWORK2 = na + IPTAU

   norma = PDLANGE('Fro',na,na,A,1,1,DESCA,WORK )
!!$   WRITE(*,*) 'The total workspace', lworkt, norma

!   CALL PDLAPRNT( na,na,A,1,1,DESCA,1,1,'A0',6, WORK )

   ttt0 = MPI_Wtime()
   call pqdwhfacs(na,na,sigma,a,1,1,desca,q,1,1,desca,WORK,lwork,descw,WORK(IPTAU),WORK(IPWORK2),&
        lwork2,info,'N' )
   ttt1 = MPI_Wtime()
   if( myid == 0) print *,'pqdwhfacs costs  :',ttt1-ttt0, 'info', info

!   CALL PDLAPRNT( na,na,Q,1,1,DESCA,1,1,'Q0',6, WORK )
   
   ttt0 = MPI_Wtime()
   CALL PQDWHSUBIT( na,A,1,1,desca,norma,Q,1,1,desca,work,lwork,work(iptau),work(ipwork2), &
        lwork2,Rk,info )
   ttt1 = MPI_Wtime()
   if( myid == 0) print *,'pqdwhsubit costs  :',ttt1-ttt0, 'info', info, 'Rk', Rk

   ! CALL PDLAPRNT( na,na,A,1,1,DESCA,1,1,'A1',6, WORK )
   !
   ! Q stores the matrix V=[V1 V2], the invariant subspaces of A.
   ! A is changed, and stores the matrix [A1 A2].
   ! 
   liwork=-1
   lwork1 = -1
   rtmp(:) =ZERO
   rrtmp(:)=0.0d0
   itmp(:) =0

   call pdsyevd('V','U',na,a,1,1,desca,work(iptau),q,1,1,desca,rtmp,lwork1,&
        itmp,liwork,info)
   liwork=itmp(1)+1

   if(myid ==0 ) print*, 'lwork, liwork ', lwork,liwork, int(rtmp)

   allocate( iwork(liwork) )

   ttt0 = MPI_Wtime()
   call pdsyevd('V','U',Rk,A,1,1,DESCA,WORK(IPTAU),WORK(IPWORK2),1,1,DESCA,WORK,lwork, &
        IWORK,liwork,info )
   ttt1 = MPI_Wtime()
   if(myid == 0) print *,'PDSYEVD costs  :',ttt1-ttt0, 'info', info, work(1)
   
   call pdgemm( 'N','N',na,Rk,Rk,ONE,Q,1,1,DESCA,WORK(IPWORK2),1,1,DESCA,ZERO,&
        WORK,1,1,DESCA )

!   If(myid ==0 ) print *, 'eigens', WORK(IPTAU:IPTAU+Rk-1)

!   CALL PDLAPRNT( na,Rk,WORK,1,1,DESCA,1,1,'QQ',6,Q )
   
   
   deallocate(a)
   deallocate(q)
   deallocate(work)
   deallocate(as)
   deallocate(iwork)

   call blacs_gridexit(my_blacs_ctxt)
!  call mpi_finalize(mpierr)
   call EXIT(STATUS)

 end program pqdwheig

!-------------------------------------------------------------------------------
